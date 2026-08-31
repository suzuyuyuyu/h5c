<h1 align="center">h5c</h1>

Serial / Parallel HDF5 を C から簡潔に扱うためのラッパーです。HDF5 C API を
直接使い、平坦バッファと明示的な shape で読み書きします。

`h5fortran` を仕様の参照元としていますが、実装上は独立しています。
C++ ラッパーは [`h5cpp`](../h5cpp) です。

## 必要なもの

- C11 対応コンパイラ
- HDF5（C コンポーネント）
- Parallel API を使う場合は Parallel HDF5 と MPI

## Build

ビルドは preset で行います。まず環境依存の設定を用意してください。

```sh
cp misc/CMakeUserPresets.json.template CMakeUserPresets.json
$EDITOR CMakeUserPresets.json      # HDF5_ROOT と CMAKE_INSTALL_PREFIX を書く
```

```sh
cmake --preset my-intel
cmake --build --preset my-intel
ctest --preset my-intel
cmake --install build/my-intel
```

`CMakePresets.json` にはコンパイラと共通オプションだけを置き、**パスなど環境に
依存する値は書きません**。それらは `CMakeUserPresets.json`（Git 管理外）で
preset を継承して足します。テンプレートが `misc/` にあります。

| preset | 用途 |
|---|---|
| `intel` / `intel-mpi` | Intel LLVM。優先して使うツールチェイン |
| `gnu` / `gnu-mpi` | 可搬性の基準。ここでも必ず通ること |

ビルドディレクトリは `build/<preset 名>` に分かれるので、`.gitignore` は
`build/` の 1 行で済みます。

Parallel は既定で OFF です。有効にすると `mpi.h` に依存する API が
`h5c/h5c_mpi.h` として現れます。Serial API は常に有効で、切り離す必要が
ないためオプションはありません。

preset を使わない場合は従来どおり指定できます。

```sh
cmake -S . -B build/manual -DHDF5_ROOT=/path/to/hdf5
```

インストール後は次の target を使います。

```cmake
find_package(h5c CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE h5c::h5c)
```

## 使い方

```c
#include <h5c/h5c.h>

double values[6] = { 1, 2, 3, 4, 5, 6 };
size_t dims[2]   = { 3, 2 };          /* 遅い次元から並べる */

h5c_file_t *f = NULL;
if (h5c_open("result.h5", H5C_TRUNCATE, &f) != H5C_OK) {
    fprintf(stderr, "%s\n", h5c_last_error()->message);
    return 1;
}
h5c_write_f64(f, "/rank/two", values, 2, dims);   /* 中間 group は自動生成 */
h5c_write_attr_str(f, "/rank/two", "units", "m/s");
h5c_close(f);
```

読み込みは **形状を問い合わせてから自分で確保する**のが基本です。

```c
h5c_dataset_info_t info;
h5c_dataset_info(f, "/rank/two", &info);   /* info.rank, dims[], type, count */

double *buf = malloc(info.count * sizeof *buf);
h5c_read_f64(f, "/rank/two", buf, info.rank, info.dims);
```

前処理・後処理では確保をライブラリに任せる経路もあります。所有権が移るので
`h5c_free()` で解放してください。

```c
double *buf = NULL;
h5c_read_alloc(f, "/rank/two", H5C_F64, (void **)&buf, &info);
/* ... */
h5c_free(buf);
```

## 配列の次元順序

**shape は row-major**、すなわち `dims[rank-1]` が最も高速に変化します。
これは HDF5 C API の既定解釈と同じで、h5c 内部での添字変換はありません。

| 言語 | 宣言 | HDF5 dims |
|---|---|---|
| Fortran | `a(nx, ny)` | `[ny, nx]` |
| C | 平坦バッファ + `dims = {ny, nx}` | `[ny, nx]` |

つまり **Fortran の `(nx, ny)` は C の `{ny, nx}`** です。`h5fortran` が書いた
ファイルはそのまま読めますし、逆も同様です。転置もコピーも発生しません。

C/C++ 側で多次元配列を持たないことを前提にしているため、rank ごとの関数は
ありません。rank は `dims` の長さで表現します。

詳細と根拠は [docs/FORMAT.md](docs/FORMAT.md) を参照してください。

## エラー処理

すべての関数が `h5c_status_t` を返します。成功は `H5C_OK`（0）です。

```c
h5c_status_t st = h5c_read_f64(f, "/maybe", buf, 1, &n);
if (st == H5C_ERR_NOT_FOUND) {
    /* 無いことが分かっているケースの処理 */
} else if (st != H5C_OK) {
    fprintf(stderr, "%s: %s\n", h5c_status_string(st),
            h5c_last_error()->message);
}
```

ファイルは **最初の非ゼロエラーを保持する sticky status** も持ちます。
戻り値を毎回見ない書き方をする場合はこちらを使ってください。

```c
h5c_status_t failed = h5c_file_status(f);   /* 何か失敗していたか */
h5c_status_t closed = h5c_close(f);         /* close は成功したか */
```

`h5c_close()` は **close 自体の成否だけ**を返します。両者を混ぜないのは、
「ファイルが flush されていない可能性がある」と「過去の処理済みの失敗」を
区別できるようにするためです。ハンドルは close で解放されるので、
sticky status は close の前に読んでください。

## 対応する型

`f32` / `f64` / `i32` / `i64` / `bool` / 文字列です。型ごとに scalar / 1D / ND の
3 本の薄いラッパーがあり、中核は `h5c_write()` / `h5c_read()` の 1 本です。

`bool` は `h5c_bool_t`（= `int8_t`）で、ファイル上は int8 基底の enum です。
1 要素 1 バイトで、`h5dump` は `TRUE` / `FALSE` と表示し、h5py は
`np.bool_` として読みます。`h5fortran` の `logical` とも往復できます。

ゼロ長 extent は正当です。要素数が 0 なら `buf` は `NULL` でよく、転送は
省略されます。Parallel でランクが 1 行も持たない場合に必要な性質で、
Serial でも同じ規則にしています。

## 多成分フィールド

ベクトル・テンソル場は `[n, ncomp]` にインターリーブして保存します
（XDMF が要求する配置です）。成分ごとに分かれた配列をそのまま渡せます。

```c
const double *comps[3] = { u, v, w };
h5c_write_interleaved_f64(f, "/fields/velocity", comps, 3, npoints);
```

書き込みはパック、読み込みは strided です。1 成分だけ欲しいときは
`h5c_read_component()` が転送量を `1/ncomp` に抑えます。理由は
`h5c/h5c.h` のコメントにあります。

## Parallel I/O

```c
#include <h5c/h5c_mpi.h>

size_t ldims[2] = { nlocal, 3 };     /* ldims[0] が分割方向。0 でもよい */

h5c_file_t *f = NULL;
h5c_popen("out.h5", H5C_TRUNCATE, &f);         /* MPI_COMM_WORLD */
h5c_pwrite(f, "/coords", local, H5C_F64, 2, ldims, H5C_WRITE_DEFAULT);

size_t offset = 0, mine = 0;
h5c_poffset(f, "/coords", &offset, &mine);     /* 自分の担当範囲 */
h5c_close(f);
```

パス `P` に書くと group ができ、`P/data` に分割方向へ連結したデータ、
`P/__partition__` に長さ `nprocs+1` のランク境界が入ります。`h5fortran` と
同一の形式です。

**すべての呼び出しは collective** で、全ランクが同じ順序で同じパスに対して
呼ぶ必要があります。引数の検証は HDF5 の collective 呼び出しより前に
`MPI_Allreduce` で全ランク集約されるため、1 ランクだけの不正な引数は
デッドロックせず全ランクで同じエラーになります。

転送は既定で collective です。`h5c_pset_collective()` で明示的にのみ
切り替わります（暗黙には変わりません）。communicator と MPI-IO ヒントを
指定する `h5c_popen_comm()` もあります。

## テスト

ラベルは「どこで実行してよいか」を表します。

```sh
ctest --preset my-intel             # 逐次テスト。ログインノードで可
sbatch scripts/run-mpi-tests.sh     # 並列テスト。バッチ投入のみ
```

test preset は `quick` のみを実行します。並列構成のビルドディレクトリで
`ctest --preset my-intel-mpi` と打っても MPI は起動しません。

**並列テストはログインノードで実行しないでください。** ランク数に関わらず
`mpiexec` の起動はバッチ経由です。そのため `mpi` ラベルは `quick` に
含めていません。ランク数はジョブスクリプトの `--rsc p=N` で変えます
（テストは任意のランク数に追従します）。

`scripts/run-mpi-tests.sh` はサイト固有の設定（partition など）を含むので、
別の環境では書き換えてください。`stdout/` と `stderr/` は Slurm が
ジョブ開始前に開くため、リポジトリに存在させています。

## Examples

`example/` に最小の使用例を置いています。**ライブラリのビルドからは参照されない
独立した CMake project** で、`find_package` で利用するため、利用者が実際に書く
コードと同じ形になっています。

```sh
cmake -S example -B example/build/intel -DCMAKE_PREFIX_PATH=...
cmake --build example/build/intel --target examples
```

詳細は [example/README.md](example/README.md) を参照してください。

## ドキュメント

- [docs/USAGE.md](docs/USAGE.md) — API の詳しい使い方
- [docs/FORMAT.md](docs/FORMAT.md) — ファイルフォーマット仕様の正本
- [docs/TODO.md](docs/TODO.md) — 今後の課題
- [log/CHANGELOG.md](log/CHANGELOG.md) — 変更履歴

各関数の契約（所有権、エラー、collective 性）は
[`include/h5c/h5c.h`](include/h5c/h5c.h) と
[`include/h5c/h5c_mpi.h`](include/h5c/h5c_mpi.h) のコメントが一次資料です。
