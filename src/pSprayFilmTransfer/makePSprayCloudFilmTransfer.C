/**
 * @file makePSprayCloudFilmTransfer.C
 * @brief 为pSpraySolutionCloud注册pSprayCloudFilmTransfer模板实例。
 * 云侧模板和膜侧fvModel共同链接在FilmTransfer库，支持运行时按字典选择。
 */
#include "pSpraySolutionCloud.H"
#include "PSprayCloudFilmTransfer.H"
makeSurfaceFilmModelType(PSprayCloudFilmTransfer, pSpraySolutionCloud);
