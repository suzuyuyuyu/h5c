# 変更履歴

## v2.0.0

4 リポジトリでバージョンを揃えた。`scheme_version` は 1 のままで、
これまでに書いた HDF5 はそのまま読める。

**破壊的変更**：`h5c_viz_open()` の引数が 5 個から 3 個になり、並列の
open は `h5c_viz_mpi.h` の `h5c_viz_popen()` に分かれた。既存の並列
呼び出しは書き換えが必要。

- serial / parallel を別ターゲットに分離し、一つのインストールで両 API を提供。可視化実装は共有し、並列 open の宣言を `h5c_viz_mpi.h` に分離。
- 通常 dataset を保存時の rank 分割に依存せず読み込む `h5c_pread_rows()` を追加。
- `test_crosslang` の h5fortran 生成物の位置を `H5C_H5FORTRAN_ARTIFACT` で指定できるようにし、未指定時のスキップを明示するようにした。従来は絶対パス固定で、他のチェックアウトでは常にスキップされていた。

## v0.1.0

- `h5fortran` を仕様の参照元とする、HDF5 C API 上の独立した C ラッパーを実装。
- ファイル open / close と明示的な3種類の mode、open mode とは独立した dataset 置換を採用。
- 平坦バッファと明示的な shape による汎用 read / write、scalar / 1D / ND / インターリーブの型別 API（f32 / f64 / i8 / i16 / i32 / i64 / bool）を追加。
- `h5c_dataset_info()`、`h5c_read_alloc()`、固定長・可変長文字列 dataset の I/O を追加。
- dataset / group / root group の文字列・数値スカラー属性を追加。
- 数値配列属性の `h5c_write_attr_array()` / `h5c_read_attr_array()` と、スカラーにも使える `h5c_attr_length()` を追加し、Fortran / C++ との相互運用を確認。
- `__partition__` による Parallel I/O、communicator / MPI-IO ヒント指定、collective（既定）/ independent の切り替えを追加。
- `[n, ncomp]` の多成分 I/O（書き込みはパック、読み込みは strided）、`h5c_poffset()` / `h5c_ppartition()` / `h5c_is_parallel()` を追加。
- 空ランクには `H5Sselect_none()` を使い、検証エラーを `MPI_Allreduce` で全ランク集約する方式を採用。
- status enum、ファイル単位の sticky status、thread-local な `h5c_last_error()` を追加し、HDF5 エラースタック出力を既定で抑制。
- bool の保存形式に int8 基底の enum を採用し、Fortran の logical との相互運用を確認。
- `h5c_viz.h` に scheme 1 の可視化 writer を追加し、h5fortran と共通レイアウトで h5xdmf による XDMF3 生成を確認。
- unstructured mesh / point cloud、複数 mesh、Scalar / Vector / Tensor6 / Tensor に対応し、型の違いは enum で扱う API を採用。
- connectivity のランクローカル 0-origin から global ID への変換と、範囲外 index の拒否を追加。
- 可視化の逐次 open を `h5c_viz_open(path, time, &out)`、並列 open を `h5c_viz_popen(path, time, comm, info, &out)` に分離（既存の並列呼び出しは移行が必要）。
- `H5C_VERSION` / `_MAJOR` / `_MINOR` / `_PATCH` と、製品バージョンから独立した `H5C_SCHEME_VERSION` を追加。
- Tensor6 の成分順を ParaView の `XX, YY, ZZ, XY, YZ, XZ` に統一。
- Fortran logical（int32）の bool 読み込みを修正し、読み込みには `H5T_NATIVE_INT8` を使用。
- 逐次インターリーブで `n == 0` の NULL 成分を許可し、存在しないファイルの open は `H5C_ERR_NOT_FOUND` を返すよう修正。
- Fortran `(nx, ny)` と C `{ny, nx}` の次元順序・バイト配置の相互運用を、転置・コピーなしで確認。
- `test_crosslang` に次元順序・型・bool enum・SPACEPAD・インターリーブ・ゼロ長 extent の検証を追加し、参照 HDF5 は実行時生成とした。
- `test_pviz` に ID 変換・空ランク・複数 mesh・connectivity の整数型・属性の検証を追加し、4ランクまで確認。
- MPI を使うテストを `mpi` ラベルに分離して `quick` から除外し、ジョブスクリプト経由で実行する運用を採用。`stdout/` / `stderr/` を配置。
