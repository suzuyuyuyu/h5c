# h5c 可視化writer実装指示

## 目的

`h5fortran`の`t_phdf5_writer`と同じHDF5 scheme 1を、h5cから出力できるようにする。
XDMF/XML生成は実装せず、独立ツール`h5xdmf`に任せる。

## 参照元

- `h5fortran/src/fypp/parallel/h5fort_parallel_visualization.fypp`
- `h5fortran/docs/USAGE-visualization.md`
- `h5xdmf/src/h5xdmf/schemes/v1.py`
- `h5xdmf/docs/design.md`

`h5fortran`のAPIを機械的に翻訳せず、ファイル上のgroup、dataset、属性、shape、dtypeを
互換にする。h5xdmf、h5fortranへのbuild依存やsubmoduleは追加しない。

## 必須条件

- `scheme_version`はh5cの製品SemVerから独立した定数`1`とする。
- root属性`scheme_version`と`time`を出力する。
- meshごとにgeometry、任意のconnectivity、point/cell dataを出力する。
- `topology_type`、`nodes_per_element`、`attribute_type`をscheme 1どおり付与する。
- Vector/Tensorは既存の`h5c_write_interleaved()`を再利用し、`[n, ncomp]`で保存する。
- `ncomp`は1、3、6、9。Tensor6は`XX, XY, XZ, YY, YZ, ZZ`順とする。
- Parallelではrank-local connectivityを0-originで受け取り、global node offsetを加える。
- connectivity範囲外、属性不整合、rank間で一致すべきmetadataの不一致は書き込み前に失敗させる。
- h5cpp向けの重複実装は作らず、h5cppからラップできる小さなC APIにする。

## 実装方針

既存のfile、attribute、parallel、interleaved APIを組み合わせ、可視化固有コードだけを
追加する。新しい依存ライブラリ、plugin機構、設定ファイル、XDMF生成機能は追加しない。

## 検証

小さなscheme 1ファイルを生成するテストを1つ追加し、HDF5 C APIでgroup、属性、shape、
値を検査する。可能ならPATH上の`h5xdmf validate`による検証を任意テストとして追加する。
MPIテストは必ず既存のジョブスクリプト経由で実行し、`mpiexec`を直接実行しない。

変更後はREADME、`docs/USAGE.md`、`docs/FORMAT.md`、`docs/TODO.md`、
`log/CHANGELOG.md`を必要な範囲だけ更新する。
