# 源码索引

各库接口、单位、源项符号和调用顺序见[库开发说明](../docs/库开发说明.md)。

| 目录 | 输出或用途 |
|---|---|
| librSolutionThermo/ | liblibrSolutionThermo.so |
| pSpraySolutionCloud/ | libpSpraySolutionCloud.so |
| pSprayFilm/ | libpSprayFilmSolver.so，含面积守恒热边界 |
| pSprayFilmPhaseChange/ | libpSprayFilmPhaseChange.so |
| pSprayFilmTransfer/ | libpSprayFilmTransfer.so |
| pSprayImplicitFilmPhaseChange/ | libpSprayImplicitFilmPhaseChange.so |
| pSpraySmoothFilmPhaseChange/ | libpSpraySmoothFilmPhaseChange.so |
| pSprayRebasedTemperature/ | libpSprayRebasedTemperature.so |
| librThermo/ | 共享物性内核与原始关联式实现 |
| pSprayInjection/ | 编入SolutionCloud的PSprayFlatFanInjection.H |

`Make/files`给出编译单元和库名，`Make/options`给出Foundation 14接口依赖。Allwmake按依赖顺序编译八库。生成的Make平台目录、lnInclude和build目录由Git忽略。
