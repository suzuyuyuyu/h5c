# Examples

用途ごとに独立した最小例を置く。**ライブラリのビルドからは参照されない**。
テストのように必須にはせず、使い方をコンパクトに残すことが目的である。

| ディレクトリ | 内容 | 実行ファイル |
|---|---|---|
| `serial/` | 読み書き、属性、形状問い合わせ、エラー処理 | `example_serial` |
| `interleaved/` | 成分ごとの配列をベクトル場として保存する | `example_serial_interleaved` |
| `parallel/` | 分散読み書きと `__partition__` | `example_parallel` |
| `parallel-interleaved/` | 分散されたベクトル場（成分別配列 + 領域分割） | `example_parallel_interleaved` |
| `visualization/` | ParaView 向け時系列。XDMF 生成まで通す | `example_parallel_visualization` |

## ビルド

`example/` は**独立した CMake project** である。`find_package(h5c)` で
インストール済みのライブラリを使うので、**利用者が実際に書く CMake と
同じ形**になっている。

```sh
cmake -S example -B example/build/intel \
      -DCMAKE_PREFIX_PATH="$HOME/.local/opt/intel/h5c-0.1.0"
cmake --build example/build/intel --target examples
```

`h5c` がインストールされていなければ、隣のソースツリーを
`add_subdirectory` して一緒にビルドする。新規 clone でもそのまま動く。

`parallel/` は `H5C_ENABLE_PARALLEL=ON` でビルドした `h5c` を指したときだけ
作られる。target 上の `H5C_HAVE_PARALLEL` を見て判断しているので、
インストール済みでもソースツリーからでも同じように働く。

```sh
cmake -S example -B example/build/intel-mpi \
      -DCMAKE_PREFIX_PATH="$HOME/.local/opt/intel/h5c-mpi-0.1.0" \
      -DCMAKE_C_COMPILER=mpiicx
cmake --build example/build/intel-mpi --target examples
```

## 実行

逐次の例はそのまま実行してよい。

```sh
cd example/build/intel && ./example_serial && ./example_serial_interleaved
```

**並列の例はログインノードで実行しないこと。** ランク数に関わらずバッチ経由で
投入する。リポジトリのルートから投入する。

```sh
sbatch example/run-parallel-example.sh
```

## 出力例

```text
$ ./example_serial
/mesh/coords: rank=2 dims={3, 2} count=6
values: 1 2 3 4 5 6
units: m
missing dataset -> not found

$ ./example_serial_interleaved
components:
  1 10 100
  2 20 200
  3 30 300
  4 40 400
as stored: 1 10 100 2 20 200 3 30 300 4 40 400
v only: 10 20 30 40

$ sbatch example/run-parallel-example.sh   # 4 ランク
rank 0: rows [0, 2)
rank 1: rows [2, 5)
rank 2: rows [5, 9)
rank 3: rows [9, 14)
global shape: {14, 3} from 4 ranks
```

`example_serial_interleaved` の `as stored` の行が、成分ごとの配列がファイル上では
`u0 v0 w0 u1 v1 w1 ...` とインターリーブされていることを示している。
これが XDMF の Vector attribute が要求する配置である。

## 分散されたベクトル場

`parallel-interleaved/` が実際のソルバーにいちばん近い。成分ごとの配列
（GPU に向いた持ち方）、MPI による領域分割、そして XDMF が要求する
`[total_n, 3]` の 1 データセットという 3 つを同時に満たす。

**ランク 0 はわざと 1 点も持たない。** 空の部分領域は正当であり、そのランクも
すべての collective 呼び出しに参加する。成分ポインタは `NULL` でよい。

```text
$ sbatch example/run-parallel-example.sh   # 4 ランク
rank 0: 0 points at rows [0, 0)
rank 1: 3 points at rows [0, 3)
rank 2: 4 points at rows [3, 7)
rank 3: 5 points at rows [7, 12)
global shape: {12, 3} from 4 ranks
as stored: 1000 2000 3000 1001 2001 3001 1002 2002 3002 2000 4000 6000 ...
```

`as stored` は、ランク 1 の 1 点目 `(u,v,w) = (1000, 2000, 3000)`、2 点目
`(1001, 2001, 3001)`、… と並んだあとにランク 2 の `(2000, 4000, 6000)` が続く。
成分が別々の配列であっても、ファイル上は点ごとにインターリーブされ、
各ランクのブロックが自分の offset に置かれていることが読み取れる。

なお属性は group ではなく `<path>/data` に付ける。group には
`__partition__` も入っているためである。

## 可視化まで通す

`visualization/` は 5 ステップの時系列を書く。四面体 mesh（圧力波と渦速度）と
粒子群（沈降）の 2 mesh を 1 ファイルに入れる。

```sh
BUILD_DIR=example/build/intel-mpi sbatch example/run-parallel-example.sh
```

出力を XDMF3 に変換する。`h5xdmf` は別プロジェクトである。

```sh
cd ../h5xdmf && uv sync
R=../h5c/example/build/intel-mpi/result
uv run h5xdmf "$R/seq*.h5" --metadata $R/metadata.h5 --outdir $R
```

```text
manifest: .../result/metadata.h5 (5 steps)
fluid:    .../result/fluid.xdmf
particles: .../result/particles.xdmf
```

`result/fluid.xdmf` と `result/particles.xdmf` を ParaView で開く。

**これは相互運用の検証も兼ねている。** `h5xdmf` は `h5fortran` 用に書かれた
ツールであり、それが `h5c` の出力をそのまま読めることが、両者のフォーマットが
同一であることの証明になっている。
