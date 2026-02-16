/**
 * @file ring_force.c
 * @brief   Gravitational force from a uniform ring of mass
 * @author  Mark Dodici <mark.a.dodici@gmail.com>
 *
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Ring Force$
 *
 * ======================= ===============================================
 * Authors                 Mark Dodici
 * Implementation Paper    n/a
 * Based on                e.g., Lass & Blitzer (1983)
 * C Example               n/a
 * Python Example          n/a
 * ======================= ===============================================
 *
 * Adds the gravitational force from a thin uniform ring of mass m at radius r
 * in the z=0 plane, centered on the origin. Uses complete elliptic integrals
 * K(k) and E(k) computed via the arithmetic-geometric mean.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * rfp_m (double)               Yes         Mass of the ring
 * rfp_r (double)               Yes         Radius of the ring
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * None
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Complete elliptic integral of the first kind K(k), |k| < 1, via AGM.
 * K(k) = pi / (2 * M(1, k')) with k' = sqrt(1 - k^2). */
static double ellint_K(double k)
{
    const double k2 = k * k;
    if (k2 >= 1.) return 0.; /* singularity; caller should avoid k^2 -> 1 */
    if (k2 <= 0.) return M_PI / 2.; /* K(0) = pi/2 */
    {
        double a = 1.;
        double b = sqrt(1. - k2);
        double tol = 2. * DBL_EPSILON * (1. + fabs(b));
        while (fabs(a - b) > tol) {
            double a_new = 0.5 * (a + b);
            double b_new = sqrt(a * b);
            a = a_new;
            b = b_new;
        }
        return M_PI / (2. * a);
    }
}

/* Complete elliptic integral of the second kind E(k), |k| < 1, via AGM.
 * Same iteration as K; E(k) = (pi/2) * (1 - sum 2^{n-1} c_n^2) / M(1,k'). */
static double ellint_E(double k)
{
    const double k2 = k * k;
    if (k2 >= 1.) return 0.; /* E(1)=1 but limit is delicate; avoid */
    if (k2 <= 0.) return M_PI / 2.; /* E(0) = pi/2 */
    {
        double a = 1.;
        double b = sqrt(1. - k2);
        double sum_c2 = 0.;
        double pow2 = 1.0; /* 2^{n-1} for n=1,2,... */
        double tol = 2. * DBL_EPSILON * (1. + fabs(b));
        while (fabs(a - b) > tol) {
            double c = 0.5 * (a - b);
            sum_c2 += pow2 * (c * c);
            pow2 *= 2.;
            double a_new = 0.5 * (a + b);
            double b_new = sqrt(a * b);
            a = a_new;
            b = b_new;
        }
        return (M_PI / 2.) * (1. - sum_c2) / a;
    }
}

void rebx_ring_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N)
{
    struct rebx_extras* const rebx = sim->extras;
    const double* rfp_m = rebx_get_param(rebx, force->ap, "rfp_m");
    const double* rfp_r = rebx_get_param(rebx, force->ap, "rfp_r");
    if (rfp_m == NULL || rfp_r == NULL) return;

    const double G = sim->G;
    const double m = *rfp_m;
    const double a = *rfp_r;
    if (a <= 0. || m <= 0.) return;

    /* Potential: Phi = -(2*G*M/pi) * K(k) / Q,  Q = sqrt((a+R)^2 + z^2),  k^2 = 4*a*R/Q^2.
     * Force F = -grad Phi. dK/dk = (E - (1-k^2)*K) / (k*(1-k^2)). */
    const double fac = (2. * G * m) / M_PI;

    for (int i = 1; i < N; i++) {
        const struct reb_particle p = particles[i];
        const double x = p.x;
        const double y = p.y;
        const double z = p.z;
        const double R = sqrt(x * x + y * y);

        if (R < 1e-30 && fabs(z) < 1e-30) continue; /* at origin, force undefined */

        const double Q2 = (a + R) * (a + R) + z * z;
        const double Q = sqrt(Q2);
        if (Q < 1e-30) continue;

        const double k2 = (4. * a * R) / Q2;
        if (k2 >= 1. - 1e-14) continue; /* on/near ring: singularity */
        const double k = sqrt(k2);

        const double Kk = ellint_K(k);
        const double Ek = ellint_E(k);
        const double one_minus_k2 = 1. - k2;
        double dK_dk;
        if (one_minus_k2 > 1e-30) {
            dK_dk = (Ek - one_minus_k2 * Kk) / (k * one_minus_k2);
        } else {
            dK_dk = 0.;
            fprintf(stderr, "reboundx ring_force: warning: 1 - k^2 too small (k^2=%.6g), dK/dk set to 0 for particle %d\n", k2, i);
        }

        /* d(k^2)/dR = 4*a*(a^2 - R^2 + z^2)/Q^4,  d(k^2)/dz = -8*a*R*z/Q^4 */
        const double dk2_dR = (4. * a) * (a * a - R * R + z * z) / (Q2 * Q2);
        const double dk2_dz = (-8. * a * R * z) / (Q2 * Q2);
        const double dk_dR = (k > 1e-30) ? (0.5 / k) * dk2_dR : 0.;
        const double dk_dz = (k > 1e-30) ? (0.5 / k) * dk2_dz : 0.;

        const double dQ_dR = (a + R) / Q;
        const double dQ_dz = z / Q;

        /* dPhi/dR = -fac * ( (dK_dk * dk_dR)/Q - K(k) * dQ_dR / Q^2 ) */
        double dPhi_dR = -fac * ((dK_dk * dk_dR) / Q - Kk * dQ_dR / Q2);
        double dPhi_dz = -fac * ((dK_dk * dk_dz) / Q - Kk * dQ_dz / Q2);

        /* F = -grad Phi: F_R = -dPhi_dR, F_z = -dPhi_dz; F_x = F_R * (x/R), F_y = F_R * (y/R) when R>0 */
        double F_R = -dPhi_dR;
        double F_z = -dPhi_dz;

        double ax = 0., ay = 0., az = F_z;
        if (R > 1e-30) {
            const double Rinv = 1. / R;
            ax = F_R * (x * Rinv);
            ay = F_R * (y * Rinv);
        }

        particles[i].ax += ax;
        particles[i].ay += ay;
        particles[i].az += az;
    }
}
