// For license details see LICENSE.

#pragma once

typedef int particle_id;
typedef int bond_id;

/** Non-owning interface for access to particles.
 */
struct ParticleReference {
    particle_id id;
    const double *pos, *vel;
};

/** Non-owning interface for access to bonds.
 */
struct BondReference {
    bond_id bid;
    particle_id pid;
    int npartners;
    const particle_id *partner_ids;
};
