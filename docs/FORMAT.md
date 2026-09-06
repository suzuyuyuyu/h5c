# HDF5 ファイル形式

h5c・h5cpp・h5fortran と h5xdmf が共有する保存形式です。
C API の呼び方は [USAGE.md](USAGE.md)、可視化の出力手順は
[USAGE-visualization.md](USAGE-visualization.md) を参照してください。

## 次元順序

h5c の `dims` は row-major で、最後の次元が最も高速に変化します。
Fortran の配列 `a(nx, ny)` に対応する C/C++ の shape は `{ny, nx}` です。

```text
Fortran: real(real64) :: a(2, 3)   ! 値 1..6 を reshape で格納
HDF5:    shape = (3, 2)
         (0,0): 1, 2
         (1,0): 3, 4
         (2,0): 5, 6
```

h5c に渡す配列の添字と shape は次のように対応します。

```text
dims = {d0, d1, d2}
flat index = (i0 * d1 + i1) * d2 + i2
```

## データ型

数値はリトルエンディアンで保存します。

| h5c の型 | ファイル上の型 |
|---|---|
| `H5C_F32` / `H5C_F64` | `H5T_IEEE_F32LE` / `H5T_IEEE_F64LE` |
| `H5C_I8` / `H5C_I16` | `H5T_STD_I8LE` / `H5T_STD_I16LE` |
| `H5C_I32` / `H5C_I64` | `H5T_STD_I32LE` / `H5T_STD_I64LE` |
| `H5C_BOOL` | `H5T_STD_I8LE` 基底の enum（`FALSE=0`, `TRUE=1`） |
| `H5C_STRING` | 固定長 `H5T_C_S1`（SPACEPAD）。可変長も選択可能 |

bool は 1 要素 1 バイトで、h5dump では `TRUE` / `FALSE`、h5py では
`np.bool_` として読めます。h5fortran の `logical` と相互に読み書きでき、
h5fortran が整数で保存した logical も h5c で読み込めます。

固定長文字列は h5fortran と相互運用できます。可変長文字列は h5c で読み書き
できますが、h5fortran では読めません。

real128 は h5c では非対応です。その dataset の型を問い合わせると
`H5C_TYPE_UNKNOWN` になります。

## 数値配列属性

数値配列属性は長さ `count` の 1 次元 dataspace に、dataset と同じ型で保存します。
h5fortran と相互運用できます。短い数値ベクトルに使い、大きなデータは dataset にします。

## Parallel の分割レイアウト

パス `P` は group になり、次の dataset を持ちます。

```text
P/data             全ランクのデータを分割軸方向に連結したもの
P/__partition__    ランク境界を表す int64 配列、長さ nprocs + 1
```

分割軸はファイル上の第 0 次元（Fortran では最終次元）です。
複数次元の同時分割は表現しません。

```text
__partition__[0]       == 0
ランク r の開始位置    == __partition__[r]
ランク r の行数        == __partition__[r+1] - __partition__[r]
__partition__[nprocs]  == data の第 0 次元の長さ
```

境界は単調非減少で、行数 0 も表現できます。
担当範囲の取得と読み込み条件は [並列 API](USAGE.md#parallel-io) を参照してください。

## 多成分フィールド

多成分の値は `[n, ncomp]` にインターリーブして保存します。
Fortran の `field(ncomp, n)` と同じ配置で、XDMF の Vector / Tensor に対応します。

| 成分数 | 可視化フィールドの `attribute_type` |
|---|---|
| 1 | `Scalar`（dataset は 1 次元） |
| 3 | `Vector` |
| 6 | `Tensor6` |
| 9 | `Tensor` |
| その他 | 属性なし |

Tensor6 は ParaView の対称テンソル規約 **`XX, YY, ZZ, XY, YZ, XZ`** の順です。
h5xdmf は成分を並べ替えません。

## 可視化レイアウト

`scheme_version = 1` の形式です。逐次・並列とも同じ配置になります。

```text
/                              attrs: scheme_version=1 (i32), time (f64)
/<mesh>/                       attrs: topology_type (固定長文字列), nodes_per_element (i32)
/<mesh>/geometry/nodes         (total_points, 3)
/<mesh>/geometry/connectivity  (total_cells, nodes_per_element)   PolyData では無し
/<mesh>/point_data/<field>     (total_points[, ncomp])
/<mesh>/cell_data/<field>      (total_cells[, ncomp])
```

各 field の `attribute_type` は [多成分フィールド](#多成分フィールド) に従います。
connectivity はファイル全体の 0-origin 節点番号です。
geometry は f32/f64、connectivity は i8/i16/i32/i64、field はこれらの数値型に対応します。

1 ファイルに複数 mesh を保存できます。h5xdmf は mesh ごとに XDMF を生成します。
`scheme_version` はライブラリの製品バージョンとは独立しています。
ファイルを読む際は root の `scheme_version` 属性で形式を確認してください。
