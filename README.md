# OFlibrary

面向 **OpenFOAM Foundation 14** 的八个喷淋、LiBr溶液、壁膜及气膜耦合扩展库。源码包含中文注释、统一代码格式、独立编译入口和所需物性表。

## 八个库

| 库 | 作用 |
|---|---|
| `liblibrSolutionThermo.so` | 温度与LiBr浓度相关的液体热物性 |
| `libpSpraySolutionCloud.so` | 携带温度和LiBr浓度的sprayCloud、平面扇形注入与逃逸账本 |
| `libpSprayFilmSolver.so` | 壁膜盐输运、内能方程及面积守恒气膜温度边界 |
| `libpSprayFilmPhaseChange.so` | 局部反应预测后冻结一个时间步的蒸发/吸收配对源 |
| `libpSprayFilmTransfer.so` | 液滴捕集、液膜回滴及质量/盐/能量/动量转移 |
| `libpSprayImplicitFilmPhaseChange.so` | 每次外迭代更新的实验性本构源与隐式修正 |
| `libpSpraySmoothFilmPhaseChange.so` | 平滑极薄膜激活的隐式修正相变源 |
| `libpSprayRebasedTemperature.so` | 保持真实温度和梯度算子不变的气膜边界重参数化 |

扇形注入的共享头文件直接编入SolutionCloud，不生成第九个注入库。三种相变实现供选择和对照，运行算例应选择一致的膜侧/气相侧配对类型。

## 编译

在安装Foundation 14的Linux或WSL2 Debian环境中：

```bash
git clone https://github.com/dusknobo/OFlibrary.git
cd OFlibrary
./Allwmake
```

`env.sh`查找`~/OpenFOAM/OpenFOAM-14/etc/bashrc`及常见系统安装路径，也支持指定环境：

```bash
export OPENFOAM14_BASHRC=/path/to/OpenFOAM-14/etc/bashrc
JOBS=2 ./Allwmake
```

八个动态库输出到`build/lib/`，物性查询工具table_audit输出到`build/bin/`。在同一终端使用库前执行：

```bash
source env.sh
```

源码分发不依赖历史预编译二进制，应在实际使用的OpenFOAM版本和ABI下编译。本仓库的编译入口只构建库和核验工具，不启动CFD求解。

## 目录

| 路径 | 内容 |
|---|---|
| `src/` | 八库、共享物性内核和扇形注入头 |
| `data/` | 水/LiBr物性表、单位与数据来源 |
| `docs/` | 接口、物理符号、数值处理、构建验证与整理记录 |
| `scripts/` | C++格式检查入口 |
| `licenses/`、`LICENSE` | 原OpenFOAM GPL与CoolProp MIT许可证 |
| `Allwmake`、`env.sh` | 编译和运行时库环境 |

运行时的`tableFile`或`solutionTable`应指向`data/librWaterTable.dat`的实际路径。X为LiBr质量分数，T为K，压力为绝对Pa。表范围0≤X≤0.65、278.15≤T≤353.15 K；浓溶液导热关联式还有独立温度适用范围。

## 说明与验证

[库开发说明](docs/库开发说明.md) · [源码索引](src/README.md) · [物性表](docs/物性表说明.md) · [验证记录](docs/验证记录.md) · [整理与上传步骤](docs/上传记录.md)

八库源自同一套已编译验证的扩展，C++运算代码保持原实现。Implicit/Smooth和参考温度重参数化属于实验性扩展；编译、格式及物性核验应与具体算例的守恒、时间步精度和物理闭合验证分别进行。

原OpenFOAM派生代码保留GPL-3.0-or-later标记及版权声明；CoolProp派生物性实现保留MIT声明，详见[NOTICE](NOTICE)。
