import pandas as pd
import plotly.graph_objects as go
from collections import OrderedDict

# =========================
# 1. 读取数据
# =========================
file_path = "桑基图.xlsx"   # 修改为你的文件路径
df = pd.read_excel(file_path)

# 清理列名
df.columns = df.columns.str.strip()

# 统一处理字符列
for col in ["Host", "Occurrence", "Phylum", "Genera"]:
    df[col] = df[col].astype(str).str.strip()

# 数值列处理
df["average_abundances"] = pd.to_numeric(df["average_abundances"], errors="coerce").fillna(0)

# 去掉无效数据
df = df[df["average_abundances"] > 0].copy()

# =========================
# 2. 定义新的维度顺序
# Host → Occurrence → Phylum → Genera
# =========================
dimensions = ["Host", "Occurrence", "Phylum", "Genera"]

# 为避免同名节点冲突，内部用带前缀ID
node_id_cols = {}
for col in dimensions:
    node_id_col = f"{col}_id"
    df[node_id_col] = col + "||" + df[col]
    node_id_cols[col] = node_id_col

# =========================
# 3. 构造节点
# =========================
node_ids = []
node_labels = []
node_levels = []

for col in dimensions:
    unique_vals = df[col].drop_duplicates().tolist()
    for val in unique_vals:
        node_ids.append(f"{col}||{val}")
        node_labels.append(val)
        node_levels.append(col)

# 去重并保持顺序
tmp = list(OrderedDict.fromkeys(zip(node_ids, node_labels, node_levels)))
node_ids = [x[0] for x in tmp]
node_labels = [x[1] for x in tmp]
node_levels = [x[2] for x in tmp]

node_dict = {nid: i for i, nid in enumerate(node_ids)}

# =========================
# 4. 节点颜色
# =========================
level_base_colors = {
    "Host": "#54A24B",        # 绿
    "Occurrence": "#E45756",  # 红
    "Phylum": "#4C78A8",      # 蓝
    "Genera": "#F58518"       # 橙
}

node_colors = [level_base_colors[level] for level in node_levels]

# =========================
# 5. 连线颜色映射
# 第一层按 Host 分类着色
# =========================
host_color_map = {
    "deer_enriched": "rgba(84,162,75,0.45)",
    "yak_enriched": "rgba(76,120,168,0.45)",
    "non_sig": "rgba(160,160,160,0.45)"
}
default_link_color = "rgba(160,160,160,0.30)"

# =========================
# 6. 构建 links
# =========================
source = []
target = []
value = []
link_colors = []
link_labels = []

for i in range(len(dimensions) - 1):
    left_dim = dimensions[i]
    right_dim = dimensions[i + 1]

    left_id = node_id_cols[left_dim]
    right_id = node_id_cols[right_dim]

    grouped = (
        df.groupby([left_dim, right_dim, left_id, right_id], as_index=False)["average_abundances"]
        .sum()
    )

    for _, row in grouped.iterrows():
        source.append(node_dict[row[left_id]])
        target.append(node_dict[row[right_id]])
        value.append(row["average_abundances"])

        # 第一层按 Host 着色，后续层级浅灰
        if left_dim == "Host":
            this_color = host_color_map.get(row[left_dim], default_link_color)
        else:
            this_color = "rgba(140,140,140,0.25)"

        link_colors.append(this_color)

        link_labels.append(
            f"{left_dim}: {row[left_dim]} → {right_dim}: {row[right_dim]}"
            f"<br>average abundances: {row['average_abundances']}"
        )

# =========================
# 7. 节点悬停信息
# =========================
node_customdata = [f"{lvl}<br>{lab}" for lvl, lab in zip(node_levels, node_labels)]

# =========================
# 8. 绘图
# =========================
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    valueformat=".0f",
    node=dict(
        pad=22,
        thickness=22,
        line=dict(color="rgba(50,50,50,0.35)", width=0.6),
        label=node_labels,
        color=node_colors,
        customdata=node_customdata,
        hovertemplate="%{customdata}<extra></extra>"
    ),
    link=dict(
        source=source,
        target=target,
        value=value,
        color=link_colors,
        customdata=link_labels,
        hovertemplate="%{customdata}<extra></extra>"
    )
)])

fig.update_layout(
    title=dict(
        text="微生物群落组成桑基图<br><sup>路径：Host → Occurrence → Phylum → Genera；连线宽度表示 average abundances</sup>",
        x=0.5,
        xanchor="center",
        font=dict(size=22, family="Microsoft YaHei, Arial", color="#222222")
    ),
    font=dict(
        size=13,
        family="Microsoft YaHei, Arial",
        color="#222222"
    ),
    width=1450,
    height=850,
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(l=30, r=30, t=90, b=30)
)

# =========================
# 9. 保存与显示
# =========================
output_html = "美化版桑基图_Host_Occurrence_Phylum_Genera.html"
fig.write_html(output_html)
fig.show()

print(f"桑基图已保存为: {output_html}")
