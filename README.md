# zsy-atf-biosensor — TetR aTF 胆酸生物传感器定向进化建模

TetR 家族变构转录因子（aTF, 204 aa）胆酸生物传感器的定向进化**计算建模**课题公开镜像。
本仓只含**进度文档、方法与可复现分析脚本、分析产物**。

## 这里有什么
- 根控制面：`PROJECT_STATE.md`（状态真相源）、`CLAIMS.md`（冻结主张台账）、`PROJECT_BACKGROUND_PROTOCOL.md`（背景+规程）、`NEXT_STAGE_ASSESSMENT_*.md`（下一阶段实验设计）。
- `dashboard/`：PI 全景汇报 + Mermaid 流程图。
- `docs/planning/`、`docs/analysis/`：方向规划、方法论、零样本/进化图分析。
- `work/`：方法与分析脚本（landscape 构建、进化图、零样本特征、序列关、推荐器）+ `work/results/` 的分析产物（`.md` / `.json`）。

## 这里没有（有意排除，不公开）
- **原始湿实验测序数据**（`*_fastq`、负筛 reads、基因型景观 csv）——未发表，归属合作方。
- 模型权重、大文件、第三方文献 PDF、组会 PPT。
- 运维脚本与任何凭据（SSH 主机/密码/IP/token）。

## 脱敏声明
公开前已做脱敏：移除真实姓名、工号、内部绝对路径与全部凭据。文中 "师兄" / "collab" 指数据提供方。

## 读序
`PROJECT_STATE.md` → `CLAIMS.md` → `PROJECT_BACKGROUND_PROTOCOL.md` → `dashboard/PI_FULL_REPORT.md` → 按需 `docs/` + `work/results/`。

## 状态一句话
结构 + 单配体的零样本机制解释已立住（Q164R 共享开关、pocket↔HTH 33Å 长程别构）；
"跨胆酸/跨蛋白泛化" 尚未建立（bootstrap CI 跨零），等配对表型（round-0 基线 + FACS）湿实验数据。
GraphWalks 任务化是当前不依赖湿实验、可推进论文的最有价值方向。
