# 浓度相关LiBr液体热物性

本库通过OpenFOAM Foundation 14热物性接口读取XLiBr质量分数，使用共同的T-X表计算密度、黏度、焓和内能。盐输运由使用该热物性的求解器负责。

```text
h(p,T,X) = hEq + (p-pEq)/rho
e(T,X) = hEq-pEq/rho
```

焓及其温度、浓度导数采用一致的Hermite插值；Cp/Cv从实际能量函数求导，The用有界Newton/二分从内能反求温度。

mixture中的tableFile指向`data/librWaterTable.dat`。X范围0—0.65，T范围278.15—353.15 K；实际浓溶液导热另要求T≥293.15 K。导热模式支持ASHRAE526RP及组件核验专用的constantForVerification。数值尾迹标记只影响有界导热参考，不修改真实T/e/rho/X。

浓度场应在初始或重启目录提供；越界物性状态直接报错。体积注入需要显式组分源，当前sourceCellComposition接口不支持默认纯水占位。

在仓库根目录执行`./Allwmake`。详细接口见[库开发说明](../../docs/库开发说明.md)，表结构见[物性表说明](../../docs/物性表说明.md)。
