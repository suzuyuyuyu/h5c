<h1 align="center">h5c</h1>

Serial / Parallel HDF5 を C から簡潔に扱うためのライブラリです。
平坦バッファと明示的な shape で読み書きし、可視化用の出力にも対応します。
C++ からは [h5cpp](https://github.com/suzuyuyuyu/h5cpp/tree/main) を利用できます。

## 必要なもの

- C11 対応コンパイラ、CMake 3.20 以上
- HDF5（C コンポーネント）
- Parallel API を使う場合は Parallel HDF5 と MPI

## Build とインストール

テンプレートの `HDF5_ROOT` と `CMAKE_INSTALL_PREFIX` を環境に合わせます。

```sh
cp misc/CMakeUserPresets.json.template CMakeUserPresets.json
$EDITOR CMakeUserPresets.json
cmake --preset my-intel
cmake --build --preset my-intel
cmake --install build/my-intel
```

GNU では `my-gnu`、並列版では `my-intel-mpi` を使います。
Parallel API は `H5C_ENABLE_PARALLEL=ON` で有効になります（既定は `OFF`）。
インストール先には `${HOME}/.local/opt/intel/h5c-<version>` などを指定してください。

利用側の CMake では、インストール先を `CMAKE_PREFIX_PATH` に指定します。

```cmake
find_package(h5c CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE h5c::h5c)
```

## 使い方

```c
#include <h5c/h5c.h>

int main(void) {
    double values[] = {1, 2, 3};
    h5c_file_t *f = NULL;
    if (h5c_open("result.h5", H5C_TRUNCATE, &f) != H5C_OK) return 1;
    h5c_status_t written = h5c_write_f64_1d(f, "/values", values, 3);
    h5c_status_t closed = h5c_close(f);
    return written != H5C_OK || closed != H5C_OK;
}
```

- [API の使い方](docs/USAGE.md)：読み込み、所有権、エラー、並列 I/O
- [可視化出力](docs/USAGE-visualization.md)：メッシュとフィールドの出力
- [ファイル形式](docs/FORMAT.md)：次元順序、型、他言語との相互運用
- [使用例](example/README.md)：用途別のプログラムとビルド方法
