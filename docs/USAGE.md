# h5c の使い方

概要とビルド方法は [README.md](../README.md)、ファイルフォーマットの仕様は
[FORMAT.md](FORMAT.md) を参照してください。この文書は API の使い方を扱います。

各関数の契約（所有権、エラー、collective 性）は `include/h5c/h5c.h` と
`include/h5c/h5c_mpi.h` のコメントが一次資料です。C ライブラリでは利用者が
必ずヘッダを読むため、そこを厚くしています。

## 初期化

`h5c_init()` は必須ではありません。どの関数も初回呼び出し時に遅延初期化します。
明示的に呼ぶ意味があるのは、HDF5 が stderr に吐く巨大なエラースタックを
抑制するタイミングを制御したい場合です。

```c
h5c_init();                        /* 任意 */
h5c_set_error_verbosity(1);        /* HDF5 のエラースタックも見たいとき */
...
h5c_finalize();                    /* 任意。H5close() を呼ぶ */
```

`h5fortran` と違い、`h5open_f` / `h5close_f` に相当する呼び出しを利用者に
義務付けていません。HDF5 C API が自動初期化されるためです。

## ファイル

mode は 3 値で、**常に明示**します。

```c
h5c_file_t *f = NULL;
h5c_status_t st = h5c_open("result.h5", H5C_TRUNCATE, &f);
```

| mode | 意味 | `h5fortran` の対応 |
|---|---|---|
| `H5C_READ` | 既存ファイルを読み取り専用で開く | `H5FORTRAN_READ_ONLY` |
| `H5C_READWRITE` | 既存ファイルを読み書きで開く | mode 省略 |
| `H5C_TRUNCATE` | 新規作成し、既存ファイルを切り詰める | `H5FORTRAN_FORCE_WRITE` |

`h5fortran` は mode 省略で read-write になりますが、これは Fortran の optional
引数があってこその設計です。C では省略という概念を持たないほうが明快なので、
3 値の enum を必ず書かせています。

既存 HDF5 コードに組み込む場合は、自分で管理している `hid_t` を借用できます。

```c
h5c_file_t *wrapper = NULL;
h5c_file_from_hid(my_fid, &wrapper);   /* 借用。close しない */
...
h5c_close(wrapper);                    /* wrapper だけ解放。my_fid は生きている */
```

逆方向は `h5c_file_hid()` です。

## 書き込み

中間 group は自動生成されます。

```c
double values[6] = { 1, 2, 3, 4, 5, 6 };
size_t dims[2]   = { 3, 2 };

h5c_write_f64(f, "/rank/two", values, 2, dims);
h5c_write_f64_1d(f, "/rank/one", values, 3);
h5c_write_f64_scalar(f, "/scalar/value", 42.5);
```

型ごとに `scalar` / `1d` / ND の 3 本があります。中核は `h5c_write()` の 1 本で、
上記はその inline ラッパーです。flags が必要なときは中核を直接呼びます。

既存 dataset への書き込みは、**形状が一致すればその場に書きます**。
形状が違う場合は `H5C_ERR_SHAPE_MISMATCH` になり、置き換えるには
`H5C_WRITE_REPLACE` を渡します。

```c
h5c_write(f, "/rank/one", other, H5C_F64, 1, &n, H5C_WRITE_REPLACE);
```

`h5fortran` は open mode の値を dataset 置換にも流用していますが、
h5c では概念を分離しています。

## 読み込み

**形状を問い合わせてから自分で確保する**のが基本形です。ソルバーでは
受け側バッファが既に確保済みであることが普通で、ライブラリが `malloc` を
握るのは速度と所有権の両面で不利なためです。

```c
h5c_dataset_info_t info;
if (h5c_dataset_info(f, "/rank/two", &info) != H5C_OK) { /* ... */ }

/* info.rank, info.dims[], info.type, info.count */
double *buf = malloc(info.count * sizeof *buf);
h5c_read_f64(f, "/rank/two", buf, info.rank, info.dims);
```

形状は**厳密に一致していなければなりません**。次元を入れ替えて渡すと
`H5C_ERR_SHAPE_MISMATCH` になります。黙って解釈し直すことはしません。

前処理・後処理では確保を任せられます。所有権が移るので `h5c_free()` で
解放してください。

```c
double *buf = NULL;
h5c_read_alloc(f, "/rank/two", H5C_F64, (void **)&buf, &info);
h5c_free(buf);
```

空の dataset でも `*buf` は非 NULL になります。`NULL` + `H5C_OK` は失敗と
区別できないためです。

存在確認は `h5c_exists()` ですが、これは「存在しない」と「ハンドルや path が
不正」の両方で 0 を返します。区別が必要なら `h5c_dataset_info()` の status を
見てください。

## 型

| `h5c_type_t` | C の型 | 備考 |
|---|---|---|
| `H5C_F32` | `float` | |
| `H5C_F64` | `double` | |
| `H5C_I32` | `int32_t` | |
| `H5C_I64` | `int64_t` | `h5fortran` の汎用 API には対応なし |
| `H5C_BOOL` | `h5c_bool_t`（`int8_t`） | ファイル上は int8 基底の enum |
| `H5C_STRING` | `char *` | 専用 API を使う |

`bool` には `h5c_bool_t` と `H5C_TRUE` / `H5C_FALSE` を使ってください。

```c
h5c_bool_t flags[4] = { H5C_TRUE, H5C_FALSE, H5C_FALSE, H5C_TRUE };
h5c_write_bool(f, "/flags", flags, 2, (size_t[]){ 2, 2 });
```

0 / 1 以外の値は検証されずそのまま書かれます。正規化には一時バッファと
全要素の走査が必要で、速度優先の方針に反するためです。

### ゼロ長 extent

要素数 0 は正当です。`buf` は `NULL` でよく、転送は省略されます。

```c
size_t empty[2] = { 0, 3 };
h5c_write_f64(f, "/empty", NULL, 2, empty);   /* {0,3} の dataset ができる */
```

Parallel でランクが 1 行も持たない場合に必要な性質で、Serial でも同じ規則に
しています。ライブラリ内に検証規則を 2 つ持たないためです。

## 文字列

```c
h5c_write_string(f, "/text/title", "flow field", H5C_WRITE_DEFAULT);

char *value = NULL;
h5c_read_string(f, "/text/title", &value);
printf("%s\n", value);
h5c_free_string(value);
```

書き込みの既定は**固定長**（`H5T_C_S1` + SPACEPAD）で、`h5fortran` が
そのまま読める形式です。読み込みは固定長・可変長の両方を受け付け、
右側の空白は取り除かれます。

可変長で書きたい場合は明示的に呼びますが、**`h5fortran` は読めません**。

```c
h5c_write_string_vlen(f, "/text/note", "...", H5C_WRITE_DEFAULT);
```

## 属性

対象は dataset、group、root group（`"/"`）です。dataset 書き込みとは
API を分けています。可変長引数や構造体配列は C では読みにくいためです。

```c
h5c_write_attr_str(f, "/rank/two", "units", "m/s");
h5c_write_attr_str(f, "/", "created_by", "solver v2");

double t = 0.125;
h5c_write_attr_scalar(f, "/", "time", &t, H5C_F64);

char *units = NULL;
h5c_read_attr_str(f, "/rank/two", "units", &units);
h5c_free_string(units);

double got = 0.0;
h5c_read_attr_scalar(f, "/", "time", &got, H5C_F64);
```

同名の属性は置き換えられます。配列属性は用途が限られるため対応していません。

## 多成分フィールド（ベクトル・テンソル）

ソルバーが `u`, `v`, `w` を別々に持っていても、そのまま渡せます。
ファイル上は `[n, ncomp]` にインターリーブされます。

```c
const double *comps[3] = { u, v, w };
h5c_write_interleaved_f64(f, "/fields/velocity", comps, 3, npoints);

double *out[3] = { gu, gv, gw };
h5c_read_interleaved_f64(f, "/fields/velocity", out, 3, npoints);
```

インターリーブが必要なのは **XDMF3 の Vector / Tensor attribute がこの配置を
要求する**からです。ParaView にベクトルとして認識させることで、
velocity magnitude の色付けやテンソル不変量の計算が使えます。

`ncomp` は 1 / 3 / 6 / 9 が `Scalar` / `Vector` / `Tensor6` / `Tensor` に対応します。
`ncomp = 6`（対称 3×3）の成分順は XDMF3 の規約に従い
**`XX, XY, XZ, YY, YZ, ZZ`** です。

1 成分だけ欲しいときは、他の成分を一切読みません。

```c
h5c_read_component_f64(f, "/fields/velocity", v_only, 1, npoints);
```

**書き込みはパック、読み込みは strided** です。stride 付きの collective 書き込みは
MPI-IO で read-modify-write を誘発するため避けています。読み込みにはその問題が
なく、転送量が `1/ncomp` で済みます。

大きなフィールドは自動でタイル分割されます。上限は変更できます。

```c
h5c_set_pack_limit(64u * 1024u * 1024u);   /* 既定は 256 MiB */
```

### 注意

`n` は**検証できません**。各成分は長さを持たない生ポインタなので、実際の長さより
大きい `n` を渡すと範囲外アクセスになります。C ではこれ以上できないので、
長さを持ち運ぶラッパー（`h5cpp` は span を使っています）を挟むのが安全です。

`ncomp == 0` はエラーですが `n == 0` は正当です。後者はそのランクが行を
持たないことを意味し、成分ポインタも `NULL` でよいという意図的な非対称です。

## エラー処理

```c
h5c_status_t st = h5c_read_f64(f, "/maybe", buf, 1, &n);
if (st == H5C_ERR_NOT_FOUND) {
    /* 無くてもよいケース */
} else if (st != H5C_OK) {
    const h5c_error_t *e = h5c_last_error();
    fprintf(stderr, "%s: %s (herr=%ld)\n",
            h5c_status_string(st), e->message, e->hdf5_err);
}
```

`h5c_last_error()` は thread-local で、次に同じスレッドで失敗するまで有効です。

ファイルは **最初の非ゼロエラーを保持する sticky status** も持ちます。
戻り値を毎回見ない書き方をする場合に使ってください。

```c
h5c_status_t failed = h5c_file_status(f);   /* 何か失敗していたか */
h5c_status_t closed = h5c_close(f);         /* close は成功したか */
h5c_file_clear_status(f);                   /* 意図的な失敗の後で消す */
```

**`h5c_close()` は close 自体の成否だけを返します。** 両者を混ぜると
「ファイルが flush されていない可能性がある」と「過去の処理済みの失敗」を
区別できず、前者に対処できなくなるためです。ハンドルは close で解放されるので、
sticky status は close の前に読んでください。

status の数値は追記のみで、既存の値は変わりません。

## Parallel I/O

`h5c/h5c_mpi.h` を include します。**このヘッダだけが `mpi.h` を含みます**。
逐次利用者は `h5c/h5c.h` だけを使い、MPI に一切依存しません。

```c
#include <h5c/h5c_mpi.h>

size_t ldims[2] = { nlocal, 3 };     /* ldims[0] が分割方向。0 でもよい */

h5c_file_t *f = NULL;
h5c_popen("out.h5", H5C_TRUNCATE, &f);                     /* MPI_COMM_WORLD */
h5c_pwrite(f, "/coords", local, H5C_F64, 2, ldims, H5C_WRITE_DEFAULT);
h5c_close(f);
```

communicator と MPI-IO ヒントを指定する形もあります。

```c
h5c_popen_comm("out.h5", H5C_READWRITE, my_comm, my_info, &f);
```

### 分割の規則

分割方向は **`dims[0]`**（最も低速に変化する軸）です。Fortran の最終次元と
同じファイル軸なので、`h5fortran` との相互運用がそのまま成立します。

`dims[0]` はランクごとに違ってよく、0 も許されます。それ以外の次元は
全ランクで一致していなければなりません。

### 自分の担当範囲を知る

```c
size_t offset = 0, mine = 0;
h5c_poffset(f, "/coords", &offset, &mine);

size_t count = 0;
h5c_ppartition(f, "/coords", NULL, 0, &count);       /* 長さを問い合わせる */
int64_t *bounds = malloc(count * sizeof *bounds);
h5c_ppartition(f, "/coords", bounds, count, NULL);   /* 全ランクの境界 */
```

`h5c_pdataset_info()` は形状（ローカルとグローバル）を返しますが、
**開始位置は返しません**。`__partition__` を自分で開く必要はありません。

### collective の規律

**すべての呼び出しは collective** で、全ランクが同じ順序で同じパスに対して
呼ぶ必要があります。

引数の検証は HDF5 の collective 呼び出しより前に `MPI_Allreduce` で
全ランク集約されます。そのため 1 ランクだけの不正な引数は**デッドロックせず、
全ランクで同じエラー**になります。これがないと、失敗したランクが早期に抜けて
残りが collective 呼び出しの中で待ち続けます。

転送は既定で collective です。切り替えは明示的にのみ行われます。

```c
h5c_pset_collective(f, 0);   /* independent へ */
```

これは通信を伴わないローカルな状態ですが、**全ランクで同じ値でなければ
なりません**。collective な転送は全ランクが同じモードを要求している必要が
あるためです。

`h5c_is_parallel()` は `h5c/h5c.h` にあるので、`mpi.h` を含まないコードからも
並列ハンドルかどうかを判定できます。

### 分散インターリーブ

```c
const double *comps[3] = { u, v, w };
h5c_pwrite_interleaved(f, "/fields/velocity", (const void *const *)comps,
                       3, nlocal, H5C_F64, H5C_WRITE_DEFAULT);
```

タイル分割が必要な場合、**タイル数は `MPI_Allreduce(MAX)` で全ランクで
合意されます**。ローカルな `n` から各ランクが独自にタイル数を決めると、
collective 呼び出しの回数が食い違ってデッドロックするためです。
行を使い切ったランクも空選択で全ての呼び出しに参加します。

## テストの実行

```sh
ctest --preset my-intel             # 逐次テスト。ログインノードで可
sbatch scripts/run-mpi-tests.sh     # 並列テスト。バッチ投入のみ
```

**並列テストはログインノードで実行しないでください。** ランク数に関わらず
`mpiexec` の起動はバッチ経由です。`mpi` ラベルを `quick` に含めていないのは、
ログインノードで習慣的に打つコマンドが誤って MPI ジョブを起動しないように
するためです。
