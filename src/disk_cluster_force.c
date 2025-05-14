/** * @file disk_cluster_force.c
 * @brief   Force resultant from a disk/cluster potential
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
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Disk+Cluster Potential$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 M. Dodici
 * Implementation Paper    n/a
 * Based on                Kaur & Stone (2025)
 * C Example               n/a
 * Python Example          n/a
 * ======================= ===============================================
 * 
 * Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
 * Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
 * and will act on all other particles.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * Acentral (double)             Yes         Normalization for central acceleration.
 * gammacentral (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_disk_cluster_force(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    const double* gc = rebx_get_param(rebx, force->ap, "dcp_gc");
    const double* gd = rebx_get_param(rebx, force->ap, "dcp_gd");
    const double* mo = rebx_get_param(rebx, force->ap, "dcp_mo");
    const double* md = rebx_get_param(rebx, force->ap, "dcp_md");
    const double* mc = rebx_get_param(rebx, force->ap, "dcp_mc");
    const double* rh = rebx_get_param(rebx, force->ap, "dcp_rh");
    const double* smooth = rebx_get_param(rebx, force->ap, "dcp_smooth");

    const double fac1 = *mo/pow(*rh,3); 
    const double tgd = (2-*gd)*(3-*gd); 
    const double C = 6*(3-*gd)/(8-tgd);
    const double B = (tgd - 2)/(2*(6-tgd));
    const double A = -(tgd-2)*(8-tgd)/(2*tgd*(6-tgd));

    for (int i=1; i<N; i++){
        const double* m = rebx_get_param(rebx, particles[i].ap, "dcp_m");
        const struct reb_particle p = particles[i];
        const double x = p.x;
        const double y = p.y;
        const double z = p.z;
        const double r2 = x*x + y*y + z*z;
        const double r = sqrt(r2);
        const double cos_t = z/r;
        const double cos_t2 = pow(cos_t,2);
        const double abs_cos_t = abs(cos_t);

        const double fac2 = *mc*pow(r/(*rh), -(*gc)); 
        const double fac3 = C * (*md/(*mo)) * pow(r/(*rh), -(*gd));
        const double fac31 = (A + abs_cos_t + B*cos_t2)*(2-(*gd));
        const double fac32xy = abs_cos_t + 2*B*cos_t2;

        const double zsmooth = 1.;
        if (smooth == 1){
            const double zsmooth = (1.- exp(-pow(z,2)/(.001*r2)));
        }

        const double sgnz = 1.;
        if (z < 0){const double sgnz = -1.;}
        
        const double fac32z = pow(sin(acos(z/r)),2) * (zsmooth*sgnz + 2*B*cos_t);

        particles[i].ax += -(fac1 * x*(fac2 + fac3*(fac31 - fac32xy)))/(*m);
        particles[i].ay += -(fac1 * y*(fac2 + fac3*(fac31 - fac32xy)))/(*m);
        particles[i].az += -(fac1 * (fac2*z + fac3*(fac31*z + fac32z*r)))/(*m);
    }
}