# h5c ファイルフォーマット仕様

このドキュメントは **C と Fortran に共通のファイルフォーマット契約**の正本である。
`h5c` と `h5fortran` は独立した実装だが、ここに書かれたレイアウトを共有することで、
一方が書いた HDF5 をもう一方が正しく読める。

Fortran API そのものの仕様は `h5fortran/docs/SPEC.md` を参照。

## 次元順序

これがもっとも間違えやすい箇所である。

**HDF5 の Fortran ライブラリは dims を反転して記録する。** `h5fortran` は
`shape(array)` をそのまま `h5screate_simple_f` に渡しているが、ファイル上の
dataspace は反転した並びになる。実測（`h5fortran` のテスト出力）で確認した結果:

```text
Fortran:  real(real64) :: a(2, 3)   ! 値 1..6 を reshape で格納
ファイル: DATASPACE SIMPLE { ( 3, 2 ) }
          (0,0): 1, 2
          (1,0): 3, 4
          (2,0): 5, 6
```

結果としてファイルは HDF5 公式が推奨する自然な形になっている。したがって **`h5c` は
転置もコピーも次元反転も行わない**。

| 言語 | 宣言 | 要素アクセス | HDF5 dims |
|---|---|---|---|
| Fortran | `a(nx, ny)` | `a(i, j)` | `[ny, nx]` |
| C / C++ | 平坦バッファ（`nx*ny` 要素） | `a[(j-1)*nx + (i-1)]` | `[ny, nx]` |

`h5c` の `dims` は **row-major**、すなわち `dims[rank-1]` が最も高速に変化する。
これは HDF5 C API の既定解釈と一致する。

```text
rank 3, dims = {d0, d1, d2}
flat index = (i0 * d1 + i1) * d2 + i2
```

column-major の実行時切り替えは提供しない。転置コピーが必須になり、
大規模並列 I/O で許容できないためである。

## データ型

| `h5c_type_t` | ファイル上の型 | メモリ上の型 | Fortran 対応 |
|---|---|---|---|
| `H5C_F32` | `H5T_IEEE_F32LE` | `float` | `real(real32)` |
| `H5C_F64` | `H5T_IEEE_F64LE` | `double` | `real(real64)` |
| `H5C_I32` | `H5T_STD_I32LE` | `int32_t` | `integer(int32)` |
| `H5C_I64` | `H5T_STD_I64LE` | `int64_t` | （h5fortran に対応なし） |
| `H5C_BOOL` | int8 基底の enum | `int8_t` | `logical` |
| `H5C_STRING` | 固定長 `H5T_C_S1` | `char*` | `character` |

ファイル上の型を明示的にリトルエンディアンにしているのは、
異なるプラットフォーム間でファイルが再現可能になるようにするためである。

### bool

ファイル上は **int8 を基底とする enum**（`FALSE=0`, `TRUE=1`）とする。

```text
DATATYPE  H5T_ENUM {
   H5T_STD_I8LE;
   "FALSE"            0;
   "TRUE"             1;
}
```

- 1 要素 1 バイト。HDF5 の要素はバイト境界に配置されるため、これが下限である
  （`H5Tset_precision()` で 1 ビット精度を宣言しても要素は 1 バイトを占める）。
- `h5dump` が `TRUE` / `FALSE` と表示する。
- h5py が `np.bool_` として読む。これは h5py が numpy の bool を保存するときの
  表現そのものである。
- `h5fortran` はこれを `logical` として読める。`h5fort_read_lgc_*` は rank しか
  検査せず `H5T_NATIVE_INTEGER` で読むため、HDF5 が自動変換する（検証済み）。

`h5fortran` は `logical` を `H5T_STD_I32LE` で書く。`h5c` はこれを `H5C_BOOL`
として読める。**ただし自動ではない。**

HDF5 は **enum → integer の変換は行うが、integer → enum の変換経路を持たない**
（実測で確認）。したがって int32 の dataset を bool の enum をメモリ型として
読もうとすると `H5Dread` が失敗する。そのため `h5c` は次の非対称な扱いをしている。

| 方向 | メモリ型 | ファイル型 |
|---|---|---|
| 書き込み | enum | enum（変換なし、memcpy） |
| 読み込み | `H5T_NATIVE_INT8` | enum でも整数でも可 |

`h5c_bool_t` は `int8_t` なので、読み込みを int8 経由にしても呼び出し側からは
何も変わらない。この非対称性は「メモリ型を 1 つに揃える」という一見無害な整理で
壊れるため、`test/test_crosslang.c` が明示的に検査している。

なお当初この文書は「HDF5 が I32→I8 変換を行うので問題ない」と書いていたが、
これは誤りだった。実際に検証していたのは逆方向（`h5fortran` が `h5c` の enum を
`H5T_NATIVE_INTEGER` で読む）であり、それを反対方向にも一般化していた。

0 / 1 以外の値は検証されずそのまま書かれる。正規化には一時バッファと
全要素の走査が必要で、速度優先の方針に反するためである。

### 文字列

書き込みの既定は固定長 `H5T_C_S1`（`h5fortran` が読める形式）。読み込みは
固定長・可変長の両方を受け付ける。可変長で書きたい場合は
`h5c_write_string_vlen()` を明示的に呼ぶが、`h5fortran` はそれを読めない。

## Parallel の分割レイアウト

path `P` へ分散配列を書くと **group** ができる。

```text
P/data             全ランクのデータを分割軸方向に連結したもの
P/__partition__    ランク境界を表す int64 配列、長さ nprocs + 1
```

```text
__partition__[0]        == 0
ランク r の開始位置     == __partition__[r]
ランク r のローカル長   == __partition__[r+1] - __partition__[r]
全体長                  == __partition__[nprocs]
```

値は単調非減少で、ローカル長 0 も表現できる。

**分割軸はファイル上の第 0 次元**である。C では `dims[0]`（最も低速に変化する軸）、
Fortran では最終次元にあたる。両者は同じファイル軸を指すため、`__partition__` は
そのまま相互運用できる。

分割軸以外の次元は全ランクで一致していなければならない。不一致は HDF5 の
collective 呼び出しより前に、全ランクで合意した上で拒否される。

`__partition__` が 1 次元の境界配列である以上、複数次元の同時分割は表現できない。

### partition へのアクセス

`__partition__` を利用者が直接開く必要はない。**レイアウトに手を伸ばさずに済む
アクセサを用意している。**

| 知りたいこと | 使う関数 |
|---|---|
| 自分の担当ブロックの開始位置と行数 | `h5c_poffset()` |
| 全ランクの境界配列 | `h5c_ppartition()` |
| ローカル・グローバルの形状 | `h5c_pdataset_info()` |

`h5c_pdataset_info()` は形状だけを返し、**開始位置は返さない**。
そのため以前は `"<path>/__partition__"` を逐次リーダーで開き、group 相対の名前を
ハードコードするしかなかった。`H5C_PARTITION_NAME` が public なのは、その時期の
名残りである。新しく書くコードでは `h5c_poffset()` を使う。

いずれも collective であり、`__partition__` の検証は `h5c_pread()` と同一である。

なお **h5cpp は現在 `__partition__` を自前で読んでいる**（`read_replicated()` 経由）。
これらのアクセサを使うように簡素化できる。

### 読み込み時の検証

読み込み前に以下をすべて検査し、いずれかが破れていればデータを転送せずに失敗する。

- `__partition__` の長さが現在の MPI プロセス数 + 1 と一致する
- 先頭要素が 0
- 単調非減少である
- 最終要素が `data` の第 0 次元の長さと一致する

## 多成分フィールド（ベクトル・テンソル）

ベクトル・テンソル場は **インターリーブして `[n, ncomp]`** で保存する。
これは XDMF3 の Vector / Tensor attribute が要求する配置であり、ParaView に
ベクトルとして認識させるために必要である（velocity magnitude の色付けや
テンソル不変量の計算がこれに依存する）。

`h5fortran` の可視化 writer は `field(ncomp, nlocal)` を受けるが、
HDF5 Fortran ライブラリの次元反転により、ファイル上は同じ `[n, ncomp]` になる。

| `ncomp` | XDMF `attribute_type` |
|---|---|
| 1 | `Scalar` |
| 3 | `Vector` |
| 6 | `Tensor6` |
| 9 | `Tensor` |

### テンソル成分の順序

`ncomp = 6`（対称 3×3）の成分順は **XDMF3 の規約に従い `XX, XY, XZ, YY, YZ, ZZ`**
とする。これに一致していなければテンソル不変量の計算が壊れる。

この順序は `h5fortran` にも `h5xdmf` にも文書化されていない
（`model.py` は成分数から型を推定するだけ）。

## 可視化レイアウト（scheme_version = 1）

`h5c_viz.h` が書く形式。`h5fortran` の `t_phdf5_writer` と同一であり、
Python の `h5xdmf` がどちらの出力からも XDMF3 を生成する。

```text
/                              attrs: scheme_version=1 (i32), time (f64)
/<mesh>/                       attrs: topology_type (固定長文字列), nodes_per_element (i32)
/<mesh>/geometry/nodes         (total_points, 3)
/<mesh>/geometry/connectivity  (total_cells, nodes_per_element)   PolyData では無し
/<mesh>/point_data/<field>     (total_points[, ncomp])
/<mesh>/cell_data/<field>      (total_cells[, ncomp])
```

`ncomp > 1` の field には `attribute_type` 属性が付く。

| `ncomp` | `attribute_type` |
|---|---|
| 1 | なし（1 次元 dataset。XDMF は Scalar と解釈する） |
| 3 | `Vector` |
| 6 | `Tensor6`（成分順 `XX, XY, XZ, YY, YZ, ZZ`） |
| 9 | `Tensor` |
| その他 | なし（XDMF に名前がないため） |

1 ファイルに複数 mesh を置ける。`h5xdmf` は mesh group ごとに `.xdmf` を出す。

**connectivity はランクローカルな 0-origin で渡し、writer が global ID に変換する。**
各ランクは自分が所有する cell だけを書く（ghost cell を書くと重複する）。
node は cell が参照する全ローカル node を含み、ランク境界での重複は許容される。
このため呼び出し側に node 番号の統合通信は不要である。

geometry は f32/f64、connectivity は i8/i16/i32/i64、field はそれらすべてに対応する。
`h5fortran` は real128 にも対応するが、**`h5c` は対応しない**。Fortran の real128 は
IEEE binary128 で、C では `long double` ではなく `__float128`（コンパイラ拡張）が
必要になるためである。real128 の dataset は `H5C_TYPE_UNKNOWN` として読まれる。

## 製品バージョンと scheme

製品バージョンは `h5c` / `h5fortran` / `h5cpp` で独立した SemVer とする。
共有するのはこのフォーマット契約だけである。

`h5fortran` は `scheme_version = 製品 major` としているが、**`h5c` は踏襲しない**。
scheme は共有する契約であり、`h5c` 自身の都合で major が上がった瞬間に
相互運用が壊れるためである。`H5C_SCHEME_VERSION`（`h5c_version.h`）は明示的に
1 とする。**reader はこの値を仮定せず、対象ファイルの `scheme_version` 属性を
読んで判断する。**

## 相互運用の検証方法

参照 `.h5` ファイルはリポジトリに含めない。代わりに、**`h5c` を一切使わず
素の HDF5 C API だけで参照ファイルを組み立てるプログラム**をテストに含め、
実行時に生成する。逆方向は参照 reader が `h5c` の出力の dataspace 次元と
バイト列を直接検査する。

こうすることで差分の読めないバイナリを版管理せずに済み、しかも期待値が
ソースコードとして明示されるため、この文書の実行可能な裏付けになる。

テストデータは**転置しても区別がつくもの**を使う。対称な形状（`2×2`）や
対称な値は、次元順序を間違えても偶然一致してしまうため使わない。
