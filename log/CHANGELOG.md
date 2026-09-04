# CHANGELOG

## 0.1.0 (未リリース)

`h5c` の最初の実装。`h5fortran` を仕様の参照元として、HDF5 C API 上に
独立した C ラッパーを構築した。

### 追加

**可視化 writer（`h5c_viz.h`、scheme_version = 1）**

- `h5fortran` の `t_phdf5_writer` と同一レイアウトの HDF5 時系列を書く。
  Python の `h5xdmf` が `h5c` の出力から XDMF3 を生成できることを実測で確認した
  （`h5fortran` 用に書かれたツールがそのまま読めることが、形式が同一である証明）。
- unstructured mesh と point cloud、1 ファイル複数 mesh、scalar / Vector /
  Tensor6 / Tensor の field。
- connectivity はランクローカル 0-origin で受け、writer が node offset を加えて
  global ID にする。範囲外 index は拒否する。呼び出し側に node 番号の統合通信は不要。
- `h5fortran` が fypp で約 60 手続きに展開している部分が 10 関数に収まる。
  型は enum が吸収するため。

**バージョン定数（`h5c_version.h`、CMake から生成）**

- `H5C_VERSION` / `_MAJOR` / `_MINOR` / `_PATCH` と `H5C_SCHEME_VERSION`。
- `h5fortran` は `scheme_version = 製品 major` としているが踏襲しない。scheme は
  共有する形式契約であり、`h5c` 自身の都合で major が上がると相互運用が壊れる。

**中核**

- ファイル操作（`h5c_open` / `h5c_close`、3 値の mode を常に明示）
- 汎用 read / write（平坦バッファ + 明示的な shape、rank は `dims` の長さで表現）
- 型ごとの薄いラッパー（scalar / 1D / ND / インターリーブ × `f32` / `f64` / `i8` / `i16` / `i32` / `i64` / `bool`）
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

### 修正

実装と検証の過程で見つかった、いずれも h5fortran との相互運用に関わるもの。

- **`h5c_read_bool()` が `h5fortran` の `logical`（int32 の 0/1）を読めなかった。**
  HDF5 は enum → integer の変換は行うが **integer → enum の変換経路を持たない**。
  bool の enum をメモリ型にも使っていたため `H5Dread` が失敗していた。読み込みだけ
  `H5T_NATIVE_INT8` を経由するようにした（書き込みは enum のまま、変換なし）。
  `docs/FORMAT.md` は「HDF5 が変換するので問題ない」と書いていたが誤りで、
  実際に検証していたのは逆方向だった。**片方向で測って両方向を主張していた。**
- **逐次のインターリーブが `n == 0` で NULL 成分を拒否していた。** ヘッダは
  「`n == 0` なら成分ポインタは NULL でよい」「逐次にも同様に当てはまる」と
  書いていたが、並列側だけ実装が追随していた。C++ では
  `std::vector<T>{}.data()` が NULL なので、空フィールドの自然な書き方が
  例外になっていた。
- **`h5c_open()` が存在しないファイルに `H5C_ERR_HDF5` を返していた。** ヘッダは
  `H5C_ERR_NOT_FOUND` を「no such file」と定義しており、呼び出し側が「無い」と
  「あるが使えない」を区別できなかった。

### 相互運用性

`h5fortran` が書いたファイルはそのまま読め、逆も同様である。
HDF5 の Fortran ライブラリが dims を反転して記録するため、
**Fortran の `(nx, ny)` は C の `{ny, nx}`** に対応し、転置もコピーも不要である。

実測で確認した（`h5fortran` の `r64_2d(2,3)` はファイル上 `{3,2}`、
h5c が同じバイト列を生成する）。

### テスト

- `test_crosslang` — 素の HDF5 C API だけで参照ファイルを組み立て、次元順序・
  ディスク上の型・bool の enum メンバ・文字列の SPACEPAD・インターリーブの
  バイト配置・ゼロ長 extent を固定する。**参照 `.h5` はコミットしない。**
  `h5fortran` の実出力があればそれとも突き合わせる。
- `test_pviz` — 可視化 writer。ランクローカル → global ID 変換、空ランク、
  複数 mesh、connectivity の全整数型、属性を検証（4 ランクで実行確認）。

### 既知の制約

- `chunking` と圧縮は未対応
- 複数次元の同時分割は `__partition__` の形式上表現できない
- 配列属性は未対応（大きなデータは dataset にすべき）
- Parallel の文字列 I/O は未対応（`h5fortran` も同様）
- 並列テストは 4 ランクまで実行して確認した
- real128 に非対応（`h5fortran` は対応）。Fortran の real128 は IEEE binary128 で、
  C では `__float128` というコンパイラ拡張が必要になるため

### 開発上の注意

- 並列テストは **ラベル `mpi`** で、`quick` には含まれない。ログインノードで
  `ctest -L quick` を打っても MPI は起動しない
- 並列テストは `sbatch scripts/run-mpi-tests.sh` で投入する
- `stdout/` と `stderr/` はリポジトリに存在させている。Slurm がジョブ開始前に
  これらを開くため、無いとジョブが 1 秒で失敗し出力も残らない
