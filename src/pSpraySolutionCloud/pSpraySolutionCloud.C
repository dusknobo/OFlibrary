/**
 * @file pSpraySolutionCloud.C
 * @brief 注册溶液parcel和cloud类型，供两侧喷雾云按同一类型创建。
 */
#include "pSpraySolutionCloud.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(pSpraySolutionParcel, 0);
addNamedToRunTimeSelectionTable(parcelCloud, pSpraySolutionCloud, thermo, pSpraySolutionCloud);
}
