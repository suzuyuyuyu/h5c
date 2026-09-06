# 可視化出力

メッシュとフィールドを HDF5 に書き、[h5xdmf](https://github.com/suzuyuyuyu/h5xdmf/blob/main/README.md) で
ParaView 用の XDMF を生成します。保存形式は [FORMAT.md](FORMAT.md#可視化レイアウト)、
リンクする target は [USAGE.md](USAGE.md#リンクするターゲットとヘッダー) を参照してください。

## 出力の手順

1 ステップにつき 1 ファイルを作り、`time` に時刻を渡します。
既存ファイルは切り詰められます。出力先ディレクトリは事前に作成してください。
以下は呼び出し順の例です。各関数の戻り値を確認してください。

```c
#include <h5c/h5c_viz.h>

h5c_viz_t *viz = NULL;
h5c_viz_mesh_t mesh = {0};
h5c_viz_open("result/seq000000.h5", t, &viz);

mesh.kind              = H5C_VIZ_UNSTRUCTURED;
mesh.name              = "fluid";
mesh.topology          = "Tetrahedron";
mesh.nodes_per_element = 4;
mesh.num_points        = npoints;
mesh.num_cells         = ncells;
h5c_viz_begin_mesh(viz, &mesh);

h5c_viz_write_nodes(viz, nodes, H5C_F64);
h5c_viz_write_connectivity(viz, conn, H5C_I32);
h5c_viz_write_point_data(viz, "Pressure", p, H5C_F64, 1);
h5c_viz_write_point_data(viz, "Velocity", v, H5C_F64, 3);
h5c_viz_write_cell_data(viz, "SubdomainID", id, H5C_I32, 1);

h5c_status_t failed = h5c_viz_status(viz);
h5c_status_t closed = h5c_viz_close(viz);
```

`h5c_viz_status()` は最初の失敗、`h5c_viz_close()` は close 自体の成否を返します。
close 後はハンドルを使えません。[エラー処理](USAGE.md#エラー処理) も参照してください。

## メッシュとバッファ

- 点数には、そのランクの cell が参照する全ローカル節点を含めます。
  ランク境界での節点の重複は許容されます。
- cell 数は所有する cell の数です。ghost cell は出力しません。
- connectivity はランク内の **0-origin** 節点番号を渡します。
  範囲外は拒否されます。1-origin や変換済みの global ID を渡さないでください。
- 座標は `num_points * 3` 要素を `x0, y0, z0, x1, y1, z1, ...` の順に渡します。
  connectivity は `num_cells * nodes_per_element` 要素です。
- field は点数または cell 数 × 成分数のバッファです。
  成分順は [FORMAT.md](FORMAT.md#多成分フィールド) に従います。

座標を x/y/z で別々に持つ場合は `h5c_viz_write_nodes_comps()`、
field には `h5c_viz_write_point_data_comps()` / `h5c_viz_write_cell_data_comps()`
を使い、各成分のポインタを配列で渡します。

点群は `H5C_VIZ_POLYDATA`、`num_cells = 0` にし、connectivity と cell data は書きません。
別のメッシュを追加するには、異なる `name` で `h5c_viz_begin_mesh()` を呼びます。
`name` 未指定時は非構造格子が `ugrid`、点群が `polydata` です。
非構造格子の topology 未指定時は `Hexahedron`、節点数 0 は 8 として扱います。

## 並列出力

MPI 初期化後、open を次に置き換えます。以降のメッシュ・書き込み・close は共通です。

```c
#include <h5c/h5c_viz_mpi.h>

h5c_viz_popen("result/seq000000.h5", t, MPI_COMM_WORLD, MPI_INFO_NULL, &viz);
```

全ランクが同じ順序で、同じメッシュ名・field 名を指定して操作します。
時刻・topology・成分数などの共通 metadata もそろえてください。
点数・cell 数はランクごとに異なってよく、0 でも同じ呼び出しに参加します。
`h5c_viz_status()` はローカルな参照です。
`h5c_viz_offsets()` で現在のメッシュの point / cell offset を取得できます。
communicator は close まで有効に保ち、MPI は close 後に終了してください。
MPI プログラムの実行はジョブスクリプトから行ってください。

## 後処理

[h5xdmf のインストール方法](https://github.com/suzuyuyuyu/h5xdmf/blob/main/README.md#インストール) に従って用意し、
出力した時系列を指定します。

```sh
h5xdmf "result/seq*.h5" --metadata result/metadata.h5 --outdir result
```

利用できるプログラムは [逐次例](../example/serial-viz/) と
[並列例](../example/parallel-viz/) を参照してください。
