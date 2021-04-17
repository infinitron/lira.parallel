# from https://github.com/sherpa/sherpa/issues/594
"""
Created on Fri Feb 21 12:23:24 2020
example:
import cstat
cstat.gof_cstat(x,y)
x: the model
y: the data
@author: longxi
"""
import numpy as np
from scipy.special import factorial


def gof_cstat(miu, obs):

    ce = np.zeros(np.size(miu))
    cv = np.zeros(np.size(miu))

    lnnm = np.log(obs / miu)
    where_are_nan = np.isnan(lnnm)
    where_are_inf = np.isinf(lnnm)
    lnnm[where_are_nan] = 0
    lnnm[where_are_inf] = 0
    co = miu - obs + obs * lnnm
    Co = 2 * sum(co)

    for i in range(np.size(miu)):
        if miu[i] == 0:
            ce[i] = -0.25 * miu[i] ** 3 + 1.38 * miu[i] ** 2
        elif 0 < miu[i] and miu[i] <= 0.5:
            ce[i] = (
                -0.25 * miu[i] ** 3 + 1.38 * miu[i] ** 2 - 2 * miu[i] * np.log(miu[i])
            )
        elif 0.5 < miu[i] and miu[i] <= 2:
            ce[i] = (
                -0.00335 * miu[i] ** 5
                + 0.04259 * miu[i] ** 4
                - 0.27331 * miu[i] ** 3
                + 1.381 * miu[i] ** 2
                - 2 * miu[i] * np.log(miu[i])
            )
        elif 2 < miu[i] and miu[i] <= 5:
            ce[i] = 1.019275 + 0.1345 * miu[i] ** (0.461 - 0.9 * np.log(miu[i]))
        elif 5 < miu[i] and miu[i] <= 10:
            ce[i] = 1.00624 + 0.604 / miu[i] ** 1.68
        else:
            ce[i] = 1 + 0.1649 / miu[i] + 0.226 / (miu[i] ** 2)

    where_are_nan = np.isnan(ce)
    where_are_inf = np.isinf(ce)
    ce[where_are_nan] = 0
    ce[where_are_inf] = 0

    k = np.arange(1, 5)
    for i in range(np.size(miu)):
        if 0 <= miu[i] and miu[i] <= 0.1:
            cv[i] = (
                4
                * sum(
                    np.exp(-miu[i])
                    * miu[i] ** k
                    / factorial(k)
                    * (miu[i] - k + k * np.log(k / miu[i])) ** 2
                )
                + 4 * (np.exp(-miu[i]) * miu[i] ** 0 / factorial(0) * (miu[i]) ** 2)
                - ce[i] ** 2
            )
        elif 0.1 < miu[i] and miu[i] <= 0.2:
            cv[i] = (
                -262 * miu[i] ** 4
                + 195 * miu[i] ** 3
                - 51.24 * miu[i] ** 2
                + 4.34 * miu[i]
                + 0.77005
            )
        elif 0.2 < miu[i] and miu[i] <= 0.3:
            cv[i] = 4.23 * miu[i] ** 2 - 2.8254 * miu[i] + 1.12522
        elif 0.3 < miu[i] and miu[i] <= 0.5:
            cv[i] = -3.7 * miu[i] ** 3 + 7.328 * miu[i] ** 2 - 3.6926 * miu[i] + 1.20641
        elif 0.5 < miu[i] and miu[i] <= 1:
            cv[i] = (
                1.28 * miu[i] ** 4
                - 5.191 * miu[i] ** 3
                + 7.666 * miu[i] ** 2
                - 3.5446 * miu[i]
                + 1.15431
            )
        elif 1 < miu[i] and miu[i] <= 2:
            cv[i] = (
                0.1125 * miu[i] ** 4
                - 0.641 * miu[i] ** 3
                + 0.859 * miu[i] ** 2
                + 1.0914 * miu[i]
                - 0.05748
            )
        elif 2 < miu[i] and miu[i] <= 3:
            cv[i] = (
                0.089 * miu[i] ** 3 - 0.872 * miu[i] ** 2 + 2.8422 * miu[i] - 0.67539
            )
        elif 3 < miu[i] and miu[i] <= 5:
            cv[i] = 2.12336 + 0.012202 * miu[i] ** (5.717 - 2.6 * np.log(miu[i]))
        elif 5 < miu[i] and miu[i] <= 10:
            cv[i] = 2.05159 + 0.331 * miu[i] ** (1.343 - np.log(miu[i]))
        else:
            cv[i] = 12 / (miu[i] ** 3) + 0.79 / (miu[i] ** 2) + 0.6747 / miu[i] + 2

    where_are_nan = np.isnan(cv)
    where_are_inf = np.isinf(cv)
    cv[where_are_nan] = 0
    cv[where_are_inf] = 0
    Ce = sum(ce)
    Cv = sum(cv)

    return [Co, Ce, Cv]


def get_cstat_gof(data, model, exposure):
    dat_x = data.x
    dat_y = data.y

    xlo = model.xlo
    xhi = model.xhi

    ymod = model.y

    miu = ymod * exposure * (xhi - xlo)
    obs = dat_y * exposure * (xhi - xlo)

    return gof_cstat(miu, obs)
