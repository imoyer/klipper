// DeltaXY kinematics stepper pulse time generation
//
// Copyright (C) 2024  Ilan E. Moyer <imoyer@mit.edu>
//
// Based on the Deltesian kinematics C module by Fabrice Gallet
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt
#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // struct stepper_kinematics
#include "trapq.h" // move_get_coord

struct deltaxy_stepper {
    struct stepper_kinematics sk;
    double arm_length_squared, driveline_position_x;
    double other_driveline_position_x;
    double build_area_depth;
    double hotend_eff_dist_squared, delta_along_arm_scaled, delta_perp_arm_scaled;
    int attached_to_hotend;
};

static double
deltaxy_stepper_calc_position(struct stepper_kinematics *sk, struct move *m, double move_time){
    // Returns the position of a DeltaXY stepper, based on an input move.

    struct deltaxy_stepper *ds = container_of(sk, struct deltaxy_stepper, sk);
    struct coord c = move_get_coord(m, move_time);

    double dx, dy;
    double y;

    if (ds->attached_to_hotend) {
        // this arm is attached to the hotend, straightforward
        dx = ds->driveline_position_x - c.x;

        // use the effective distance of the hotend, taking into account offset
        dy = sqrt(ds->hotend_eff_dist_squared - dx * dx);

        y = (c.y + dy) - ds->build_area_depth;
    } else {
        double dx_other = ds->other_driveline_position_x - c.x;

        // use the effective distance of the hotend, taking into account offset
        double dy_other = sqrt(ds->hotend_eff_dist_squared - dx_other * dx_other);

        // vec along arm: [-dx_other, -dy_other]
        // vec perpendicular to arm: [dy_other, -dx_other]
        double x_joint = c.x - dx_other * ds->delta_along_arm_scaled + dy_other * ds->delta_perp_arm_scaled;
        double y_joint = c.y - dy_other * ds->delta_along_arm_scaled - dx_other * ds->delta_perp_arm_scaled;

        dx = x_joint - ds->driveline_position_x;

        dy = sqrt(ds->arm_length_squared - dx * dx);

        y = (y_joint + dy) - ds->build_area_depth;
    }

    return y;
}

struct stepper_kinematics * __visible
deltaxy_stepper_alloc(double arm_length_squared, double driveline_position_x, double other_driveline_position_x, double build_area_depth, double hotend_eff_dist_squared, double delta_along_arm_scaled, double delta_perp_arm_scaled, int attached_to_hotend) {
    struct deltaxy_stepper *ds = malloc(sizeof(*ds));
    memset(ds, 0, sizeof(*ds));
    ds->arm_length_squared = arm_length_squared;
    ds->driveline_position_x = driveline_position_x;
    ds->other_driveline_position_x = other_driveline_position_x;
    ds->build_area_depth = build_area_depth;
    ds->hotend_eff_dist_squared = hotend_eff_dist_squared;
    ds->delta_along_arm_scaled = delta_along_arm_scaled;
    ds->delta_perp_arm_scaled = delta_perp_arm_scaled;
    ds->attached_to_hotend = attached_to_hotend;
    ds->sk.calc_position_cb = deltaxy_stepper_calc_position;
    ds->sk.active_flags = AF_X | AF_Y;
    return &ds->sk;
}
