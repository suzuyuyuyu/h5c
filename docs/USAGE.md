# h5c の使い方

概要とビルド方法は [README.md](../README.md)、ファイルフォーマットの仕様は
[FORMAT.md](FORMAT.md) を参照してください。この文書は API の使い方を扱います。
以下の断片例ではエラー分岐を省略しています。実際のコードでは
[エラー処理](#エラー処理) に従って戻り値を確認してください。

## リンクするターゲットとヘッダー

| ターゲット | 提供する API |
|---|---|
| `h5c::h5c_serial` | `h5c.h` の直列 I/O と `h5c_viz.h` の可視化操作 |
| `h5c::h5c_parallel` | `h5c_mpi.h` の並列 I/O と `h5c_viz_mpi.h` の `h5c_viz_popen()`。直列ターゲットに依存し、並列構成で提供 |
| `h5c::h5c` | 互換用。直列構成では serial、並列構成では serial + parallel |

```cmake
find_package(h5c CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE h5c::h5c_serial)
```

並列構成の serial ターゲットも HDF5 経由で MPI に依存します。
リンクまで MPI 不要にするには Serial HDF5 と `H5C_ENABLE_PARALLEL=OFF` を使います。
HDF5 の型 `hid_t` は公開ヘッダーから利用できます。

## 初期化

`h5c_init()` は省略できます。HDF5 のエラースタックも表示したい場合は
`h5c_set_error_verbosity(1)` を呼びます。

```c
h5c_init();                        /* 任意 */
h5c_set_error_verbosity(1);        /* HDF5 のエラースタックも見たいとき */
...
h5c_finalize();                    /* 任意。H5close() を呼ぶ */
```

## ファイル

mode は 3 値で、**常に明示**します。

```c
h5c_file_t *f = NULL;
h5c_status_t st = h5c_open("result.h5", H5C_TRUNCATE, &f);
```

| mode | 意味 |
|---|---|
| `H5C_READ` | 既存ファイルを読み取り専用で開く |
| `H5C_READWRITE` | 既存ファイルを読み書きで開く |
| `H5C_TRUNCATE` | 新規作成し、既存ファイルを切り詰める |

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

型別の `scalar` / `1d` / ND 関数を使えます。flags を指定する場合は
`h5c_write()` を使います。

既存 dataset への書き込みは、**形状が一致すればその場に書きます**。
形状が違う場合は `H5C_ERR_SHAPE_MISMATCH` になり、置き換えるには
`H5C_WRITE_REPLACE` を渡します。

```c
h5c_write(f, "/rank/one", other, H5C_F64, 1, &n, H5C_WRITE_REPLACE);
```

## 読み込み

形状を問い合わせ、必要な要素数のバッファを確保します。

```c
h5c_dataset_info_t info;
if (h5c_dataset_info(f, "/rank/two", &info) != H5C_OK) { /* ... */ }

/* info.rank, info.dims[], info.type, info.count */
double *buf = malloc(info.count * sizeof *buf);
h5c_read_f64(f, "/rank/two", buf, info.rank, info.dims);
free(buf);
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

空の dataset でも `*buf` は非 NULL になります。

存在確認は `h5c_exists()` ですが、これは「存在しない」と「ハンドルや path が
不正」の両方で 0 を返します。区別が必要なら `h5c_dataset_info()` の status を
見てください。

## 型

数値型は `float`、`double`、`int8_t`、`int16_t`、`int32_t`、`int64_t` です。
型指定には対応する `H5C_F32` / `H5C_F64` / `H5C_I8` / `H5C_I16` /
`H5C_I32` / `H5C_I64` を使います。保存形式は [FORMAT.md](FORMAT.md#データ型)、
配列の shape は [次元順序](FORMAT.md#次元順序) を参照してください。

`bool` には `h5c_bool_t` と `H5C_TRUE` / `H5C_FALSE` を使ってください。

```c
h5c_bool_t flags[4] = { H5C_TRUE, H5C_FALSE, H5C_FALSE, H5C_TRUE };
h5c_write_bool(f, "/flags", flags, 2, (size_t[]){ 2, 2 });
```

0 / 1 以外の値は検証・正規化されません。

### ゼロ長 extent

要素数 0 は正当です。`buf` は `NULL` でよく、転送は省略されます。

```c
size_t empty[2] = { 0, 3 };
h5c_write_f64(f, "/empty", NULL, 2, empty);   /* {0,3} の dataset ができる */
```

## 文字列

```c
h5c_write_string(f, "/text/title", "flow field", H5C_WRITE_DEFAULT);

char *value = NULL;
h5c_read_string(f, "/text/title", &value);
printf("%s\n", value);
h5c_free_string(value);
```

読み込みは末尾の空白を除去します。可変長で書く場合は次を使います。
保存形式と相互運用の制約は [FORMAT.md](FORMAT.md#データ型) を参照してください。

```c
h5c_write_string_vlen(f, "/text/note", "...", H5C_WRITE_DEFAULT);
```

## 属性

対象は dataset、group、root group（`"/"`）です。

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

同名の属性は置き換えられます。数値配列属性には `h5c_write_attr_array()` /
`h5c_read_attr_array()`、要素数の問い合わせには `h5c_attr_length()` を使います。
スカラー属性にも使え、要素数は1です。保存形式と用途は
[FORMAT.md](FORMAT.md#数値配列属性) を参照してください。

## 多成分フィールド（ベクトル・テンソル）

ソルバーが `u`, `v`, `w` を別々に持っていても、そのまま渡せます。
成分の配置と順序は [FORMAT.md](FORMAT.md#多成分フィールド) を参照してください。

```c
const double *comps[3] = { u, v, w };
h5c_write_interleaved_f64(f, "/fields/velocity", comps, 3, npoints);

double *out[3] = { gu, gv, gw };
h5c_read_interleaved_f64(f, "/fields/velocity", out, 3, npoints);
```

1 成分だけ欲しいときは、他の成分を一切読みません。

```c
h5c_read_component_f64(f, "/fields/velocity", v_only, 1, npoints);
```

成分配列の書き込みに使う一時バッファの上限を変更できます。

```c
h5c_set_pack_limit(64u * 1024u * 1024u);   /* 既定は 256 MiB */
```

### 注意

各成分のバッファは少なくとも `n` 要素必要です。長さは検証されません。
`ncomp == 0` はエラーです。`n == 0` なら成分ポインタは `NULL` でも構いません。

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
```

`h5c_close()` は close 自体の成否だけを返し、ハンドルを解放します。
過去のエラーは close 前に確認してください。処理済みのエラーを消すには、
ハンドルが有効な間に `h5c_file_clear_status(f)` を呼びます。

## Parallel I/O

`h5c/h5c_mpi.h` を include します。MPI の宣言はこのヘッダーと
可視化の `h5c/h5c_viz_mpi.h` に分離されています。

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

保存レイアウトは [FORMAT.md](FORMAT.md#parallel-の分割レイアウト) を参照してください。

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
free(bounds);
```

`h5c_pdataset_info()` は形状（ローカルとグローバル）を返しますが、
**開始位置は返しません**。`__partition__` を自分で開く必要はありません。

### collective の規律

並列 open・close・読み書き・形状と担当範囲の問い合わせは collective です。
全ランクが同じ順序で同じパスに対して呼びます。MPI は open 前に初期化し、
close 後に終了してください。MPI プログラムはジョブスクリプトから実行します。

転送は既定で collective です。切り替えは明示的にのみ行われます。

```c
h5c_pset_collective(f, 0);   /* independent へ */
```

`h5c_pset_collective()` はローカルな設定ですが、全ランクで同じ値にしてください。

`h5c_is_parallel()` は `h5c/h5c.h` にあるので、`mpi.h` を含まないコードからも
並列ハンドルかどうかを判定できます。

### 分散インターリーブ

```c
const double *comps[3] = { u, v, w };
h5c_pwrite_interleaved(f, "/fields/velocity", (const void *const *)comps,
                       3, nlocal, H5C_F64, H5C_WRITE_DEFAULT);
```

行数 0 のランクも同じ呼び出しに参加します。

### 読み込み条件

`h5c_pread()` は保存時と同じ MPI プロセス数で使います。
保存された境界配列の長さ・先頭・単調性・全体長が不正な場合は読み込みに失敗します。
通常の dataset から任意の連続行を読むには `h5c_pread_rows()` を使います。
これは保存済み partition を参照せず、offset とローカル shape を指定します。

## 可視化と制約

メッシュの出力手順は [USAGE-visualization.md](USAGE-visualization.md) を参照してください。
chunking・圧縮と分散文字列配列の I/O は未対応です。
