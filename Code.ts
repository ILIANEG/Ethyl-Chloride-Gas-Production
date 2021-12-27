/**
 * @customfunction
 * @param {number[][]} flowrates
 * @param {number[][]} KVALUES
 * @param {number[][]} initialVolume
 * @param {number} finalVolume
 * @param {number} molarData
 * @param {number} initialTemp
 * @param {number} alpha
 * @param {number[][]} termalProps
 * @param {number[][]} heatOfReaction
 * @param {number} reactionTemp
 * @param {number} U
 * @param {number} a
 * @param {number} T_a
 * @param {number} Error
 */
function pfr(
  // Row of molar flow rates
  F_ABCDEI_INIT: number[][],
  p: number,
  KVALUES: number[][],
  PBED: number[][],
  V_INIT: number,
  V_FINAL: number,
  MOLAR: number[][],
  T_INIT: number,
  P_INIT: number,
  Cp: number[][],
  H_R: number[][],
  T_R: number,
  U: number,
  T_a: number,
  ERR: number
): number[][] {
  let C_T0 = calcCIdeal(T_INIT, P_INIT);
  let RHO_B = PBED[3][0];
  let mu = calcViscousityAir(T_INIT);
  let rhoGas = calcRhoIdeal(F_ABCDEI_INIT, MOLAR, C_T0);
  let q = calcQFlow(F_ABCDEI_INIT, MOLAR, rhoGas);
  let ALPHA = calcAlpha(q, P_INIT, rhoGas, PBED[1][0], PBED[3][0], mu, PBED[4][0], PBED[2][0], PBED[0][0]);
  let Ua = (U * 4) / PBED[2][0];
  let [v_ii, h_ii, F_ii, T_ii, p_ii] = [V_INIT, 0.0001, F_ABCDEI_INIT[0].slice(), T_INIT, p];
  let [v_i, F_i, T_i, p_i, h_i] = [v_ii, F_ii, T_ii, p_ii, h_ii];

  while (v_ii < V_FINAL) {
    [v_i, T_i, F_i, p_i, h_i] = [v_ii, T_ii, F_ii.slice(), p_ii, h_ii];

    let k: number[][] = [[]];
    let offset: number[] = [];

    k[0] = KCalc(
      F_ABCDEI_INIT,
      F_i,
      KVALUES,
      RHO_B,
      T_INIT,
      T_i,
      p_i,
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );
    for (let i = 0; i < 8; i++) {
      offset[i] = (1 / 4) * k[0][i];
    }

    k[1] = KCalc(
      F_ABCDEI_INIT,
      F_i.map((fi, i) => fi + offset[i]),
      KVALUES,
      RHO_B,
      T_INIT,
      T_i + offset[7],
      p_i + offset[6],
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );
    for (let i = 0; i < 8; i++) {
      offset[i] = (3 / 32) * k[0][i] + (9 / 32) * k[1][i];
    }

    k[2] = KCalc(
      F_ABCDEI_INIT,
      F_i.map((fi, i) => fi + offset[i]),
      KVALUES,
      RHO_B,
      T_INIT,
      T_i + offset[7],
      p_i + offset[6],
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );
    for (let i = 0; i < 8; i++) {
      offset[i] = (1932 / 2197) * k[0][i] - (7200 / 2197) * k[1][i] + (7296 / 2197) * k[2][i];
    }

    k[3] = KCalc(
      F_ABCDEI_INIT,
      F_i.map((fi, i) => fi + offset[i]),
      KVALUES,
      RHO_B,
      T_INIT,
      T_i + offset[7],
      p_i + offset[6],
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );
    for (let i = 0; i < 8; i++) {
      offset[i] = (439 / 216) * k[0][i] - 8 * k[1][i] + (3680 / 513) * k[2][i] - (845 / 4104) * k[3][i];
    }

    k[4] = KCalc(
      F_ABCDEI_INIT,
      F_i.map((fi, i) => fi + offset[i]),
      KVALUES,
      RHO_B,
      T_INIT,
      T_i + offset[7],
      p_i + offset[6],
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );
    for (let i = 0; i < 8; i++) {
      offset[i] =
        (-8 / 27) * k[0][i] + 2 * k[1][i] - (3544 / 2565) * k[2][i] + (1859 / 4104) * k[3][i] - (11 / 40) * k[4][i];
    }

    k[5] = KCalc(
      F_ABCDEI_INIT,
      F_i.map((fi, i) => fi + offset[i]),
      KVALUES,
      RHO_B,
      T_INIT,
      T_i + offset[7],
      p_i + offset[6],
      h_i,
      ALPHA,
      C_T0,
      Cp,
      H_R,
      T_R,
      Ua,
      T_a,
      v_i,
      V_INIT
    );

    let d45: any = {
      F: {
        4: F_i.map(
          (fi, i) => fi + ((25 / 216) * k[0][i] + (1408 / 2565) * k[2][i] + (2197 / 4101) * k[3][i] - (1 / 5) * k[4][i])
        ),
        5: F_i.map(
          (fi, i) =>
            fi +
            ((16 / 135) * k[0][i] +
              (6656 / 12825) * k[2][i] +
              (28561 / 56430) * k[3][i] -
              (9 / 50) * k[4][i] +
              (2 / 55) * k[5][i])
        )
      },
      p: {
        4: p_i + ((25 / 216) * k[0][6] + (1408 / 2565) * k[2][6] + (2197 / 4101) * k[3][6] - (1 / 5) * k[4][6]),
        5:
          p_i +
          ((16 / 135) * k[0][6] +
            (6656 / 128025) * k[2][6] +
            (28561 / 56430) * k[3][6] -
            (9 / 50) * k[4][6] +
            (2 / 55) * k[5][6])
      },
      T: {
        4: T_i + ((25 / 216) * k[0][7] + (1408 / 2565) * k[2][7] + (2197 / 4101) * k[3][7] - (1 / 5) * k[4][7]),
        5:
          T_i +
          ((16 / 135) * k[0][7] +
            (6656 / 128025) * k[2][7] +
            (28561 / 56430) * k[3][7] -
            (9 / 50) * k[4][7] +
            (2 / 55) * k[5][7])
      }
    };

    let flat: number[] = [];
    let key: keyof typeof d45;
    for (key in d45) {
      if (key === "F") {
        d45[key]["5"].forEach((_: number, i: number) => {
          flat.push((1 / h_i) * Math.abs(d45[key]["5"][i] - d45[key]["4"][i]));
        });
      } else {
        flat.push((1 / h_i) * Math.abs(d45[key]["5"] - d45[key]["4"]));
      }
    }

    let maxR = Math.max(...flat);

    if (maxR < ERR) {
      F_ii = d45["F"][5];
      p_ii = d45["p"][5];
      T_ii = d45["T"][5];
      v_ii = v_i + h_i;
    }

    h_ii = (h_i * 0.84 * (ERR / maxR)) ^ (1 / 4);

    if (h_ii < 0.000001) {
      F_ii = d45["F"][5];
      p_ii = d45["p"][5];
      T_ii = d45["T"][5];
      v_ii = v_i + h_i;
      h_ii = 0.000001;
    } else if (0.0001 < h_ii) {
      F_ii = d45["F"][5];
      p_ii = d45["p"][5];
      T_ii = d45["T"][5];
      v_ii = v_i + h_i;
      h_ii = 0.0001;
    }
  }

  let solution: number[] = [];

  for (let i = 0; i < F_i.length; i++) {
    if (v_ii == V_FINAL) {
      solution.push(F_ii[i]);
    } else {
      solution.push((F_i[i] * (v_ii - V_FINAL) + F_ii[i] * (V_FINAL - v_i)) / (v_ii - v_i));
    }
  }
  if (v_ii == V_FINAL) {
    solution.push(p_ii, T_ii);
  } else {
    solution.push(
      (p_i * (v_ii - V_FINAL) + p_ii * (V_FINAL - v_i)) / (v_ii - v_i),
      (T_i * (v_ii - V_FINAL) + T_ii * (V_FINAL - v_i)) / (v_ii - v_i)
    );
  }

  return [[...solution, V_FINAL]];
}

function KCalc(
  F_ABCDEI_INIT: number[][],
  F_ABCDEI: number[],
  KVALUES: number[][],
  RHO_B: number,
  T_INIT: number,
  T: number,
  p: number,
  h_i: number,
  ALPHA: number,
  C_T0: number,
  CpData: number[][],
  Hr: number[][],
  T_r: number,
  Ua: number,
  T_a: number,
  v_i: number,
  v_init: number
): number[] {
  let fs: number[] = [];
  let r_i: number[];
  let c: number[] = [];
  //let cps: number[];

  for (let i = 0; i < F_ABCDEI.length; i++) {
    c[i] = (((C_T0 * F_ABCDEI[i]) / totalF(F_ABCDEI)) * p * T) / T_INIT;
  }
  let rr = rates(c, KVALUES, T);
  r_i = specieReactionRate(rr);
  for (let i = 0; i < F_ABCDEI.length; i++) {
    fs[i] = RHO_B * h_i * r_i[i];
  }

  fs.push((((((-h_i * ALPHA) / (2 * p)) * totalF(F_ABCDEI)) / totalF(F_ABCDEI_INIT[0])) * T) / T_INIT);
  fs.push(h_i * calcT(F_ABCDEI, Hr, CpData, rr, Ua, T, T_r, T_a, v_i, v_init));
  return fs;
}

function rates(C: number[], KVALUES: number[][], T: number) {
  let ki: number[] = [];
  for (let i = 0; i < KVALUES.length; i++) {
    ki[i] = calcKi(KVALUES[i][0], T);
  }

  let keqs = calcKeqs(T);

  return [
    ki[0] * (C[0] * C[1] - (C[2] * C[3]) / keqs[0]),
    ki[1] * (C[0] ^ (2 - (C[4] * C[3]) / keqs[1])),
    ki[2] * (C[4] * C[1] - (C[2] * C[0]) / keqs[2])
  ];
}

function specieReactionRate(r: number[]) {
  return [-r[0] - 2 * r[1] + r[2], -r[0] - r[2], r[0] + r[2], r[0] + r[1], r[1] + r[2], 0];
}

function calcT(
  Fi: number[],
  Hr: number[][],
  CpData: number[][],
  ri: number[],
  Ua: number,
  T: number,
  Tr: number,
  Ta: number,
  vi: number,
  vInit: number
) {
  let intCps = calcIntCps(CpData, T, Tr);
  let cPs = calcCps(CpData, T);
  let sumDHr =
    -(-ri[0] * (Hr[2][0] + Hr[3][0] - Hr[1][0] - Hr[0][0] + (intCps[2] + intCps[3] - intCps[1] - intCps[0]))) +
    (-2 *
      ri[0] *
      ((1 / 2) * Hr[2][0] + (1 / 2) * Hr[4][0] - Hr[0][0] + ((1 / 2) * intCps[1] + (1 / 2) * intCps[4] - intCps[0])) +
      ri[2] * (Hr[0][0] + Hr[3][0] - Hr[1][0] - Hr[4][0] + (intCps[0] + intCps[3] - intCps[1] - intCps[4])));
  let fCpSum = 0;
  Fi.forEach((f, i) => {
    fCpSum += f * cPs[i];
  });
  return (sumDHr + Ua * (Ta - T)) / fCpSum;
}

function calcKi(k: number, T: number, T_MEAN = 300) {
  return k * Math.exp(1 / T - 1 / T_MEAN);
}

function calcKeqs(T: number) {
  return [
    0.0271 * T ** 2 - 17.999 * T + 3150.2,
    0.0005 * T ** 2 - 0.32 * T + 59.375,
    0.0002 * T ** 2 - 0.2074 * T + 77.991
  ];
}

function calcCps(CpData: number[][], T: number) {
  let cps: number[] = [];
  CpData.forEach((c) => {
    cps.push((c[0] + c[1] * T + c[2] * T ** 2) * 4.184);
  });
  return cps;
}

function calcIntCps(CpData: number[][], Ti: number, Tf: number) {
  let icps: number[] = [];
  CpData.forEach((c) => {
    icps.push(((c[2] * (Tf ** 3 - Ti ** 3)) / 3 + c[1] * (Tf - Ti)) * 4.184);
  });
  return icps;
}

function totalF(F_ABCDEI: number[]) {
  return F_ABCDEI.reduce((x, y) => x + y);
}

function calcCIdeal(T: number, P: number, R = 8.314) {
  return P / (R * (T + 273.15));
}

function calcRhoIdeal(F_ABCDEI: number[][], molar: number[][], C: number) {
  let avgMolar = 0;
  F_ABCDEI[0].forEach((f: number, i: number, fs: number[]) => {
    avgMolar += (f / totalF(fs)) * molar[i][0];
  });
  return avgMolar * C;
}

function calcViscousityAir(T: number) {
  return (
    ((0.01827 * (0.555 * 524.07 + 120)) / (0.555 * (((T + 273.15) * 9) / 5) + 120)) *
    (((T + 273.15) * 9) / 5 / 524.07) ** (3 / 2) *
    0.001
  );
}

function calcQFlow(F_ABCDEI: number[][], molar: number[][], rho: number) {
  let q = 0;
  F_ABCDEI[0].forEach((f, i, fs) => {
    q += f * molar[i][0];
  });
  return q / rho;
}

function calcAlpha(
  qFlow: number,
  P: number,
  rhoGas: number,
  rhoCat: number,
  rhoBulk: number,
  mu: number,
  voidage: number,
  bedDiameter: number,
  particleDiameter: number
) {
  let uSup = qFlow / ((Math.PI / 4) * bedDiameter ** 2);
  let beta =
    ((rhoGas * uSup * (1 - voidage)) / (rhoGas * particleDiameter * voidage ** 3)) *
    ((150 * (1 - voidage) * mu) / particleDiameter + 1.75 * rhoGas * uSup);
  return (2 * beta) / ((((1 - voidage) * Math.PI) / 4) * bedDiameter ** 2 * rhoCat * P);
}
