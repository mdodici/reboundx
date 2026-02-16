/**
 * @file cluster_force.c
 * @brief   Force from the gradient of a cluster potential (power-law mass profile)
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
 * $Cluster Force$
 *
 * ======================= ===============================================
 * Authors                 M. Dodici
 * Implementation Paper    n/a
 * Based on                Cluster potential with mass-density power law
 * C Example               n/a
 * Python Example          n/a
 * ======================= ===============================================
 *
 * Adds the force from the gradient of a cluster potential: F = nabla Phi_c.
 * With D = 1/(2 - gamma_c),
 *   nabla Phi_c = (G * m_c / r_0^3) * (r / r_0)^(-gamma_c) * x,
 * where x is the position vector from the cluster center (origin) and r = |x|.
 * gamma_c is the mass-density power law exponent. gamma_c = 2 is invalid (D undefined).
 *
 * **Effect Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * cf_mc (double)               Yes         Mass scale of the cluster (m_c)
 * cf_r0 (double)               Yes         Characteristic radius (r_0)
 * cf_gamma_c (double)          Yes         Mass-density power law exponent; must not equal 2
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

void rebx_cluster_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N)
{
    struct rebx_extras* const rebx = sim->extras;
    const double* cf_mc   = rebx_get_param(rebx, force->ap, "cf_mc");
    const double* cf_r0   = rebx_get_param(rebx, force->ap, "cf_r0");
    const double* cf_gamma_c = rebx_get_param(rebx, force->ap, "cf_gamma_c");

    if (cf_mc == NULL || cf_r0 == NULL || cf_gamma_c == NULL) return;

    const double G = sim->G;
    const double m_c = *cf_mc;
    const double r_0 = *cf_r0;
    const double gamma_c = *cf_gamma_c;

    if (fabs(gamma_c - 2.) < 1e-15) {
        fprintf(stderr, "reboundx cluster_force: error: gamma_c must not equal 2 (D = 1/(2 - gamma_c) is undefined).\n");
        exit(1);
    }
    if (r_0 <= 0.) return;

    /* F = nabla Phi_c = (G * m_c / r_0^3) * (r / r_0)^(-gamma_c) * x = G * m_c * r_0^(gamma_c - 3) * r^(-gamma_c) * x */
    const double fac = G * m_c * pow(r_0, gamma_c - 3.);

    for (int i = 1; i < N; i++) {
        const struct reb_particle p = particles[i];
        const double x = p.x;
        const double y = p.y;
        const double z = p.z;
        const double r2 = x * x + y * y + z * z;

        if (r2 < 1e-30) continue;

        const double r_inv_gamma = pow(r2, -gamma_c * 0.5);
        const double ax = fac * r_inv_gamma * x;
        const double ay = fac * r_inv_gamma * y;
        const double az = fac * r_inv_gamma * z;

        particles[i].ax += ax;
        particles[i].ay += ay;
        particles[i].az += az;
    }
}
