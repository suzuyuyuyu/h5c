# CHANGELOG

## 0.1.0 (未リリース)

`h5c` の最初の実装。`h5fortran` を仕様の参照元として、HDF5 C API 上に
独立した C ラッパーを構築した。

### 追加

**中核**

- ファイル操作（`h5c_open` / `h5c_close`、3 値の mode を常に明示）
- 汎用 read / write（平坦バッファ + 明示的な shape、rank は `dims` の長さで表現）
- 型ごとの薄いラッパー（scalar / 1D / ND の 3 本 × `f32` / `f64` / `i32` / `i64` / `bool`）
- 形状問い合わせ `h5c_dataset_info()` と、確保も行う `h5c_read_alloc()`
- 文字列 dataset（固定長が既定、可変長は明示 API、読みは両対応）
- 文字列・数値スカラー属性（dataset / group / root group が対象）
- Parallel I/O（`__partition__` による分割、collective / independent、
  communicator と MPI-IO ヒントの指定）
- 多成分フィールドのインターリーブ I/O（`[n, ncomp]`、書きはパック、読みは strided）
- レイアウトアクセサ `h5c_poffset()` / `h5c_ppartition()`
- `h5c_is_parallel()`（`mpi.h` を含まないコードから分岐できる）

**エラー処理**

- 全関数が `h5c_status_t` を返す。原因分類が enum の値だけで判別できる粒度
- ファイルごとの sticky status（最初の非ゼロを保持）
- thread-local な `h5c_last_error()`（status、生の `herr_t`、メッセージ）
- HDF5 の巨大なエラースタック出力を既定で抑制

### `h5fortran` から意図的に変えた点

- **communicator を選べる**。`h5fortran` は `MPI_COMM_WORLD` 決め打ち
- **collective / independent を明示的に切り替えられる**。既定は collective
- **ゼロサイズ選択に `H5Sselect_none()` を使う**。長さ 0 の hyperslab 選択は
  HDF5 や MPI-IO の実装で扱いが揺れうるため
- **検証エラーを `MPI_Allreduce` で全ランク集約する**。1 ランクだけが早期に
  抜けて残りが collective 呼び出しへ進むデッドロックを防ぐ
- **`bool` をファイル上で int8 基底の enum にした**。`h5fortran` は int32 で、
  1 要素 4 バイトを使っている。int8 なら 1 バイトで済み、`h5dump` は
  `TRUE` / `FALSE` と表示し、h5py は `np.bool_` として読む。
  `h5fortran` はこれを `logical` として読めることを実測で確認した
- **open mode と dataset 置換を分離した**。`h5fortran` は `mode` を両方に
  流用しているが、これは概念の重ね掛けである
- **`integer(int64)` に対応した**。`h5fortran` の汎用 API にはない

### 相互運用性

`h5fortran` が書いたファイルはそのまま読め、逆も同様である。
HDF5 の Fortran ライブラリが dims を反転して記録するため、
**Fortran の `(nx, ny)` は C の `{ny, nx}`** に対応し、転置もコピーも不要である。

実測で確認した（`h5fortran` の `r64_2d(2,3)` はファイル上 `{3,2}`、
h5c が同じバイト列を生成する）。

### 既知の制約

- `chunking` と圧縮は未対応
- 複数次元の同時分割は `__partition__` の形式上表現できない
- 配列属性は未対応（大きなデータは dataset にすべき）
- Parallel の文字列 I/O は未対応（`h5fortran` も同様）
- 並列テストは 4 ランクまで実行して確認した

### 開発上の注意

- 並列テストは **ラベル `mpi`** で、`quick` には含まれない。ログインノードで
  `ctest -L quick` を打っても MPI は起動しない
- 並列テストは `sbatch scripts/run-mpi-tests.sh` で投入する
- `stdout/` と `stderr/` はリポジトリに存在させている。Slurm がジョブ開始前に
  これらを開くため、無いとジョブが 1 秒で失敗し出力も残らない
