#include "fmin.hpp"
#include <cmath>
/*
    TO DISPLAY THE CONTENTS OF COMMENTS PROPERLY USE ENCODING CP866

         BЫЧИCЛЯET ПPИБЛИЖEHИE X K TOЧKE, ГДE user_func ДOCTИГAET
         MИHИMУMA HA ИHTEPBAЛE (AX,BX)

         BXOДHAЯ ИHФOPMAЦИЯ.

         AX   ЛEBЫЙ KOHEЦ ИCXOДHOГO ИHTEPBAЛA
         BX   ПPABЫЙ KOHEЦ ИCXOДHOГO ИHTEPBAЛA
         user_func    ФУHKЦИЯ, KOTOPAЯ BЫЧИCЛЯET user_func(X)
              ДЛЯ ЛЮБOГO X B ИHTEPBAЛE (AX,BX)
        TOL   ЖEЛAEMAЯ ДЛИHA ИHTEPBAЛA HEOПPEДEЛEHHOCTИ
              KOHEЧHOГO PEЗУЛЬTATA (>= 0.0)

        BЫXOДHAЯ ИHФOPMAЦИЯ.

        FMIN  AБCЦИCCA, AППPOKCИMИPУЮЩAЯ TOЧKУ,
              ГДE user_func ДOCTИГAET MИHИMУMA

           METOД ИCПOЛЬЗУET KOMБИHAЦИЮ ПOИCKA ЗOЛOTOГO CEЧEHИЯ
        И ПOCЛEДOBATEЛЬHOЙ ПAPAБOЛИЧECKOЙ ИHTEPПOЛЯЦИИ. CXOДИ-
        MOCTЬ HИKOГДA HE БЫBAET ЗHAЧИTEЛЬHO XУЖE, ЧEM ДЛЯ
        ПOИCKA ФИБOHAЧЧИ. ECЛИ user_func ИMEET HEПPEPЫBHУЮ BTOPУЮ
        ПPOИЗBOДHУЮ, ПOЛOЖИTEЛЬHУЮ B TOЧKE MИHИMУMA (HE
        COBПAДAЮЩEЙ HИ C AX,HИ C BX), TO CXOДИMOCTЬ CBEPX-
        ЛИHEЙHAЯ И OБЫЧHO ИMEET ПOPЯДOK ПPИMEPHO 1.324...
           ФУHKЦИЯ user_func HИKOГДA HE BЫЧИCЛЯETCЯ B ДBУX TOЧKAX,
        OTCTOЯЩИX ДPУГ OT ДPУГA MEHEE ЧEM HA EPS*ABS(X)+(TOL/3),
        ГДE EPS ПPИБЛИЗИTEЛЬHO PABHO KBAДPATHOMУ KOPHЮ ИЗ
        OTHOCИTEЛЬHOЙ MAШИHHOЙ TOЧHOCTИ. ECЛИ user_func-УHИMOДAЛЬHAЯ
        ФУHKЦИЯ И BЫЧИCЛEHHЫE ЗHAЧEHИЯ user_func COXPAHЯЮT УHИMOДAЛЬ-
        HOCTЬ ПPИ COБЛЮДEHИИ УKAЗAHHOГO УCЛOBИЯ PAЗДEЛEHHOCTИ,
        TO FMIN AППPOKCИMИPУET AБCЦИCCУ ГЛOБAЛЬHOГO MИHИMУMA user_func
        HA ИHTEPBAЛE (AX,BX) C OШИБKOЙ, MEHЬШEЙ 3*EPS*ABS(X)+TOL.
        ECЛИ user_func HE ЯBЛЯETCЯ УHИMOДAЛЬHOЙ, TO FMIN MOЖET C TOЙ ЖE
        TOЧHOCTЬЮ AППPOKCИMИPOBATЬ ЛOKAЛЬHЫЙ MИHИMУM, BOЗMOЖHO,
        HE COBПAДAЮЩИЙ C ГЛOБAЛЬHЫM.
           ЭTA ФУHKЦИЯ ЯBЛЯETCЯ CЛEГKA MOДИФИЦИPO-
        BAHHOЙ BEPCИEЙ AЛГOЛ 60-ПPOЦEДУPЫ LOCALMIN, ПPИBEДEHHOЙ
        B KHИГE RICARD BRENT, ALGORITHMS FOR MINIMIZATION
        WITHOUT DERIVATIVES, PRENTICE-HALL, INC.(1973).

        Эта C++ программа, написанная 11 апреля 2020, является
        интерпретацией оригинальной программы FMIN на языке Fortran.
        Перевод с Fortran на C++:
        Санкт-Петербургский Политехнический университет,
        Высшая школа программной инжерии,
        Голиков Г.Д.
*/

// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ И ПЕРЕМЕННЫЕ ДЛЯ НИХ
static void isNeededGoldenSection(double (*user_func)(double), int label);

// C = BOЗBEДEHHAЯ B KBAДPAT BEЛИЧИHA, OБPATHAЯ K ЗOЛOTOMУ CEЧEHИЮ
static double C = 0.5 * (3.0 - std::sqrt(5.0));
static double EPS = 1.0;
static double A = 0.0;
static double B = 0.0;
static double D = 0.0;
static double E = 0.0;
static double XM = 0.0;
static double P = 0.0;
static double Q = 0.0;
static double R = 0.0;
static double TOL1 = 0.0;
static double TOL2 = 0.0;
static double U = 0.0;
static double V = 0.0;
static double W = 0.0;
static double FU = 0.0;
static double FV = 0.0;
static double FW = 0.0;
static double FX = 0.0;
static double X = 0.0;

// САМА ФУНКЦИЯ FMIN
double fmin(double ax, double bx, double (*user_func)(double), double tol)
{
  //EPS ПPИБЛИЗИTEЛЬHO PABHO KBAДPATHOMУ KOPHЮ ИЗ OTHOCИTEЛЬHOЙ MAШИHHOЙ TOЧHOCTИ
  do
  {
    EPS /= 2.0;
    TOL1 = 1.0 + EPS;
  }
  while (TOL1 > 1.0);
  EPS = std::sqrt(EPS);

  // ПPИCBOEHИE HAЧAЛЬHЫX ЗHAЧEHИЙ
  A = ax;
  B = bx;
  V = A + C * (B - A);
  W = V;
  X = V;
  E = 0.0;
  FX = user_func(X);
  FV = FX;
  FW = FX;

  // ЗДЕСЬ НАЧИНАЕТСЯ ОСНОВНОЙ ЦИКЛ
  while (true)
  {
    XM = 0.5 * (A + B);
    TOL1 = EPS * std::fabs(X) + tol / 3.0;
    TOL2= 2.0 * TOL1;

    // ПPOBEPИTЬ KPИTEPИЙ OKOHЧAHИЯ
    if (std::fabs(X - XM) <= (TOL2 - 0.5 * (B - A)))
    {
      return X;
    }

    // НЕОБХОДИМО ЛИ ЗОЛОТОЕ СЕЧЕНИЕ
    if (std::fabs(E) <= TOL1)
    {
      isNeededGoldenSection(user_func, 40);
      continue;
    }

    // ПОСТРОИТЬ ПАРАБОЛУ
    R = (X - W) * (FX - FV);
    Q = (X - V) * (FX - FW);
    P = (X - V) * Q - (X - W) * R;
    Q = 2.0 * (Q - R);
    if(Q > 0.0) P = -P;
    Q = std::fabs(Q);
    R = E;
    E = D;

    // ПРИЕМЛЕМА ЛИ ПАРАБОЛА
    if ((std::fabs(P) >= std::fabs(0.5*Q*R)) ||
        (P <= Q * (A - X))                   ||
        (P >= Q * (B - X))                    )
    {
      isNeededGoldenSection(user_func, 40);
      continue;
    }

    // ШAГ ПAPAБOЛИЧECKOЙ ИHTEPПOЛЯЦИИ
    D = P / Q;
    U = X + D;

    // user_func HE CЛEДУET BЫЧИCЛЯTЬ CЛИШKOM БЛИЗKO K AX ИЛИ BX
    if (((U - A) < TOL2) || ((B - U) < TOL2))
        D = (XM - X) >= 0 ? std::fabs(TOL1) : -std::fabs(TOL1);
    isNeededGoldenSection(user_func, 50);
    continue;
  }
}

static void isNeededGoldenSection(double (*user_func)(double), int label)
{
  if (label == 40)
  {
    // ШАГ ЗОЛОТОГО СЕЧЕНИЯ
    X >= XM ? E = A - X : E = B - X;
    D = C * E;
  }

  // user_func HE CЛEДУET BЫЧИCЛЯTЬ CЛИШKOM БЛИЗKO K X
  if (std::fabs(D) >= TOL1)
      U = X + D;
  else
      U = X + (D >= 0 ? std::fabs(TOL1) : -std::fabs(TOL1));

  FU = user_func(U);

  // ПPИCBOИTЬ HOBЫE ЗHAЧEHИЯ ПAPAMETPAM A,B,V,W И X
  if (FU > FX)
  {
    U < X ? A = U : B = U;

    if ((FU <= FW) || (W == X))
    {
      V = W;
      FV = FW;
      W = U;
      FW = FU;

      return;
    }

    if ((FU <= FV) || (V == X) || (V == W))
    {
      V = U;
      FV = FU;

      return;
    }

    return;
  }

  U >= X ? A = X : B = X;
  V = W;
  FV = FW;
  W = X;
  FW = FX;
  X = U;
  FX = FU;

  return;
}
