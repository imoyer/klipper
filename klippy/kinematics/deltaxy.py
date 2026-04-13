# Code for handling the kinematics of DeltaXY positioning systems
#
# Copyright (C) 2024  Ilan E. Moyer <imoyer@mit.edu>
#
# Based on the "delta" kinematics module by Kevin O'Connor and the "deltesian" kinematics module by Fabrice Gallet.
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging
import math

import mathutil
import stepper


class DeltaXYKinematics:
    """Adds Klipper support for the DeltaXY positioning system."""

    def __init__(self, toolhead, config):
        """Initializes the DeltaXY kinematics object.

        toolhead -- a toolhead.ToolHead object.
        config -- a configfile.ConfigWrapper object.
        """

        # READ CONFIG FILE
        left_stepper_config = config.getsection("stepper_left")
        right_stepper_config = config.getsection("stepper_right")
        z_stepper_config = config.getsection("stepper_z")
        deltaxy_config = config.getsection("deltaxy")

        # Read Arm Lengths
        self.left_arm_length_mm = left_stepper_config.getfloat(
            "arm_length_mm", above=0.0
        )  # enforce > 0. Note that arm lengths are not necessarily the same, for asymmetric mechanisms.
        self.right_arm_length_mm = right_stepper_config.getfloat(
            "arm_length_mm", above=0.0
        )

        # Read Mechanism Parameters
        self.driveline_separation_mm = deltaxy_config.getfloat(
            "driveline_separation_mm", above=0.0
        )  # the distance between the lines of motion of both shoulders
        self.build_area_width_mm = deltaxy_config.getfloat(
            "build_area_width_mm", above=0.0
        )  # horizontal width (in the x direction) of the build area
        self.build_area_depth_mm = deltaxy_config.getfloat(
            "build_area_depth_mm", above=0.0
        )  # vertical depth (in the y direction) of the build area
        self.build_area_centerline_offset_mm = deltaxy_config.getfloat(
            "build_area_centerline_offset_mm"
        )  # offset between the centerline of the build area and the midline between the drivelines

        # Read Stepper Endstop Positions
        self.left_endstop_position_mm = left_stepper_config.getfloat(
            "position_endstop", above=0.0
        )
        self.right_endstop_position_mm = right_stepper_config.getfloat(
            "position_endstop", above=0.0
        )
        self.z_endstop_position_mm = z_stepper_config.getfloat(
            "position_endstop", above=0.0
        )

        # One of the interesting properties of DeltaXY is that it permits travel outside the standard rectangular "rational" build area.
        # While we could rely solely on the travel limits of the rails, this exposes the mechanism to singularities and correspondingly large changes
        # in the transmission ratio between the motor and toolhead motions.
        # By simply clipping the actual travel envelope in the X direction, we can avoid these singularities.
        self.toolhead_limit_min_x_mm = deltaxy_config.getfloat(
            "toolhead_limit_min_x_mm"
        )
        self.toolhead_limit_max_x_mm = deltaxy_config.getfloat(
            "toolhead_limit_max_x_mm"
        )

        # Read offset of the nozzle. This is relative to the joint, in the global X and Y coordinates, when the machine is at rest (same position on both rails)
        self.toolhead_offset_x_mm = deltaxy_config.getfloat("toolhead_offset_x_mm")
        self.toolhead_offset_y_mm = deltaxy_config.getfloat("toolhead_offset_y_mm")

        # Read Maximum Z Velocity
        max_velocity_xy, max_accel_xy = toolhead.get_max_velocity()
        self.max_velocity_z = config.getfloat(
            "max_z_velocity", max_velocity_xy, above=0.0, maxval=max_velocity_xy
        )
        self.max_accel_z = config.getfloat(
            "max_z_accel", max_accel_xy, above=0.0, maxval=max_accel_xy
        )

        # INSTANTIATE ALL RAILS
        # Each "rail" corresponds to the linear motion resulting from the rotation of a stepper motor
        # Note that we use the min and max position parameters to limit travel along the rail, so we require this of the config file.
        # The left and right rails may have different endstop positions, whether due to mfg error or especially if the mechanism is assymetric.

        self.rail_left = stepper.LookupMultiRail(
            left_stepper_config, need_position_minmax=True
        )
        self.rail_right = stepper.LookupMultiRail(
            right_stepper_config, need_position_minmax=True
        )
        self.rail_z = stepper.LookupMultiRail(
            z_stepper_config, need_position_minmax=True
        )
        self.rails = [self.rail_left, self.rail_right, self.rail_z]

        # Read Rail Travel Extents
        self.left_rail_min_position_mm, self.left_rail_max_position_mm = (
            self.rail_left.get_range()
        )
        self.right_rail_min_position_mm, self.right_rail_max_position_mm = (
            self.rail_right.get_range()
        )
        self.z_rail_min_position_mm, self.z_rail_max_position_mm = (
            self.rail_z.get_range()
        )

        # SET UP ITERATIVE SOLVER
        # These are c helper functions. We feed them with as many pre-computed variables as possible to minimize their compute time.

        # Pre-Compute Values
        self.left_arm_length_squared_mm2 = self.left_arm_length_mm**2
        self.right_arm_length_squared_mm2 = self.right_arm_length_mm**2
        self.left_driveline_x_pos_mm = (
            self.build_area_width_mm - self.driveline_separation_mm
        ) / 2 - self.build_area_centerline_offset_mm
        self.right_driveline_x_pos_mm = (
            self.build_area_width_mm + self.driveline_separation_mm
        ) / 2 - self.build_area_centerline_offset_mm

        # Pre-Compute constants for offset compensation (assume hotend attached to right arm)

        angle_right_arm_to_horizontal_rest_rad = math.acos(
            (
                self.right_arm_length_squared_mm2
                + self.driveline_separation_mm**2
                - self.left_arm_length_squared_mm2
            )
            / (2 * self.right_arm_length_squared_mm2 * self.driveline_separation_mm)
        )

        delta_x_rest = self.right_arm_length_mm * math.cos(
            angle_right_arm_to_horizontal_rest_rad
        )
        delta_y_rest = self.right_arm_length_mm * math.sin(
            angle_right_arm_to_horizontal_rest_rad
        )

        delta_x_hotend = delta_x_rest - self.toolhead_offset_x_mm
        delta_y_hotend = delta_y_rest - self.toolhead_offset_y_mm

        # effective length when taking account the offset (distance from hotend to shoulder)
        self.right_arm_length_effective_mm = math.sqrt(
            delta_x_hotend**2 + delta_y_hotend**2
        )
        self.right_arm_length_effective_squared_mm2 = (
            self.right_arm_length_effective_mm**2
        )

        # coordinate system of the effective arm
        vec_shoulder_to_hotend_x = -delta_x_hotend / self.right_arm_length_effective_mm
        vec_shoulder_to_hotend_y = -delta_y_hotend / self.right_arm_length_effective_mm
        vec_shoulder_to_hotend_perp_x = (
            delta_y_hotend / self.right_arm_length_effective_mm
        )
        vec_shoulder_to_hotend_perp_y = (
            -delta_x_hotend / self.right_arm_length_effective_mm
        )

        # offsets in that coordinate system
        self.delta_hotend_along_effective_mm = vec_shoulder_to_hotend_x * (
            -self.toolhead_offset_x_mm
        ) + vec_shoulder_to_hotend_y * (-self.toolhead_offset_y_mm)
        self.delta_hotend_perpendicular_to_effective_mm = (
            vec_shoulder_to_hotend_perp_x * (-self.toolhead_offset_x_mm)
            + vec_shoulder_to_hotend_perp_y * (-self.toolhead_offset_y_mm)
        )

        # To avoid division in the chelper
        self.delta_hotend_along_effective_scaled = (
            self.delta_hotend_along_effective_mm / self.right_arm_length_effective_mm
        )
        self.delta_hotend_perpendicular_to_effective_scaled = (
            self.delta_hotend_perpendicular_to_effective_mm
            / self.right_arm_length_effective_mm
        )

        # Set Up Rail Iterative Solvers
        self.rail_left.setup_itersolve(
            "deltaxy_stepper_alloc",
            self.left_arm_length_squared_mm2,
            self.left_driveline_x_pos_mm,
            self.right_driveline_x_pos_mm,
            self.build_area_depth_mm,
            self.right_arm_length_effective_squared_mm2,
            self.delta_hotend_along_effective_scaled,
            self.delta_hotend_perpendicular_to_effective_scaled,
            0,
        )
        self.rail_right.setup_itersolve(
            "deltaxy_stepper_alloc",
            self.right_arm_length_squared_mm2,
            self.right_driveline_x_pos_mm,
            self.left_driveline_x_pos_mm,
            self.build_area_depth_mm,
            self.right_arm_length_effective_squared_mm2,
            self.delta_hotend_along_effective_scaled,
            self.delta_hotend_perpendicular_to_effective_scaled,
            1,
        )
        self.rail_z.setup_itersolve("cartesian_stepper_alloc", b"z")

        # Configure Trapezoidal Move Queue
        for this_stepper in self.get_steppers():
            this_stepper.set_trapq(toolhead.get_trapq())
            toolhead.register_step_generator(this_stepper.generate_steps)

        # REGISTER EVENT HANDLERS
        config.get_printer().register_event_handler(
            "stepper_enable:motor_off", self._motor_off
        )

        # INITIALIZE STATE FLAGS
        self.homed_axes = [False] * 3  # tracks which axes have been homed

        # INITIALIZE TOOLHEAD POSITION
        self.set_position([0.0, 0.0, 0.0], ())

        # CALCULATE TOOL HOME POSITION
        # This is defined as the tool position when the axes are homed against the endstops
        self.tool_home_position = self._actuator_to_cartesian(
            [
                self.left_endstop_position_mm,
                self.right_endstop_position_mm,
                self.z_endstop_position_mm,
            ]
        )

    def get_steppers(self):
        """Returns a list of all stepper objects."""
        return [
            this_stepper for rail in self.rails for this_stepper in rail.get_steppers()
        ]

    def _actuator_to_cartesian(self, stepper_positions_mm):
        """Returns the cartesian coordinates of the toolhead based on the positions of the actuators.

        stepper_positions_mm -- a list of stepper positions, in the form [left_stepper_position_mm, right_stepper_position_mm, z_stepper_position_mm]
        """
        left_shoulder_y_position_mm, right_shoulder_y_position_mm, z_position_mm = (
            tuple(stepper_positions_mm)
        )

        # first, find the difference in height between the two shoulders
        delta_y_mm = (
            right_shoulder_y_position_mm - left_shoulder_y_position_mm
        )  # vertical distance between shoulders

        # this and the driveline separation gives us the linear distance between the shoulders
        diagonal_shoulder_dist_mm = math.sqrt(
            self.driveline_separation_mm**2 + delta_y_mm**2
        )

        # Next, we find the angle between the left-going horizontal and the diagonal line between shoulders
        # note that we are treating angles below the horizontal as positive for this calculation
        angle_diagonal_to_horizontal_rad = math.atan(
            delta_y_mm / self.driveline_separation_mm
        )

        # we then use the law of cosines to find the angle between the right arm and the shoulder diagonal
        angle_right_arm_to_diagonal_rad = math.acos(
            (
                self.right_arm_length_squared_mm2
                + diagonal_shoulder_dist_mm**2
                - self.left_arm_length_squared_mm2
            )
            / (2 * self.right_arm_length_mm * diagonal_shoulder_dist_mm)
        )

        # Now we finally get what we're after, which is the angle between the right arm and vertical
        angle_right_arm_to_vertical_rad = (math.pi / 2) - (
            angle_diagonal_to_horizontal_rad + angle_right_arm_to_diagonal_rad
        )

        # From here it's trivial to find the cartesian position of the toolhead relative to the right shoulder
        right_shoulder_to_tool_x_distance_mm = self.right_arm_length_mm * math.sin(
            angle_right_arm_to_vertical_rad
        )
        right_shoulder_to_tool_y_distance_mm = self.right_arm_length_mm * math.cos(
            angle_right_arm_to_vertical_rad
        )

        # Last, we transpose this into build area coordinates
        tool_position_x_mm = (
            self.right_driveline_x_pos_mm - right_shoulder_to_tool_x_distance_mm
        )
        tool_position_y_mm = (
            right_shoulder_y_position_mm
            + self.build_area_depth_mm
            - right_shoulder_to_tool_y_distance_mm
        )

        # To take offset into account, compare the rotation to a resting pose
        angle_right_arm_to_horizontal_rest_rad = math.acos(
            (
                self.right_arm_length_squared_mm2
                + self.driveline_separation_mm**2
                - self.left_arm_length_squared_mm2
            )
            / (2 * self.right_arm_length_mm * self.driveline_separation_mm)
        )

        angle_right_arm_to_vertical_rest_rad = (
            math.pi / 2
        ) - angle_right_arm_to_horizontal_rest_rad

        delta_angle = (
            angle_right_arm_to_vertical_rest_rad - angle_right_arm_to_vertical_rad
        )

        # cordinate system of the offset after rotation from the resting pose
        vec_dy_x = -math.sin(delta_angle)
        vec_dy_y = math.cos(delta_angle)
        vec_dx_x = math.cos(delta_angle)
        vec_dx_y = math.sin(delta_angle)

        # apply X offset
        tool_position_x_mm += self.toolhead_offset_x_mm * vec_dx_x
        tool_position_y_mm += self.toolhead_offset_x_mm * vec_dx_y

        # apply Y offset
        tool_position_x_mm += self.toolhead_offset_y_mm * vec_dy_x
        tool_position_y_mm += self.toolhead_offset_y_mm * vec_dy_y

        return [tool_position_x_mm, tool_position_y_mm, z_position_mm]

    def calc_position(self, stepper_position_dict):
        """Returns the cartesian coordinates of the toolhead based on the positions of the actuators.

        stepper_position_dict -- a dictionary in the form {rail_name: stepper_position_mm, ...} Note: I'm assuming this, but could be wrong about the format
        Returns the toolhead position in the format [x_mm, y_mm, z_mm]
        """
        stepper_position_list_mm = [
            stepper_position_dict[rail.get_name()] for rail in self.rails
        ]
        return self._actuator_to_cartesian(stepper_position_list_mm)

    def set_position(self, newpos, homing_axes):
        """I don't fully understand the mechanics of this function, so just using the template format.

        homing_axes -- a list of axis indices that have been homed.
        """
        for rail in self.rails:
            rail.set_position(newpos)

        for axis_index in homing_axes:
            self.homed_axes[axis_index] = True

    def home(self, homing_state):
        """Homes the axes of the machine."""

        homing_axes = homing_state.get_axes()
        home_xy = (0 in homing_axes) or (
            1 in homing_axes
        )  # we home in XY if either X or Y axes request homing
        home_z = 2 in homing_axes

        if home_xy:
            homing_state.set_axes([0, 1])
            home_positions = self.tool_home_position[:2] + [None, None]
            force_positions = [None, -self.build_area_depth_mm * 0.5, None, None]
            homing_state.home_rails(self.rails[:2], force_positions, home_positions)

        if home_z:
            homing_state.set_axes([2])
            home_positions = [None, None, self.tool_home_position[2], None]
            force_positions = [None, None, -self.tool_home_position[2], None]
            homing_state.home_rails([self.rails[2]], force_positions, home_positions)

    def _motor_off(self, print_time):
        self.homed_axes = [False] * 3

    def check_move(self, move):
        """Raises a move error if the move is outside of the bounds of motion for the machine.

        Unlike linear positioning mechanisms, DeltaXY has a weirdly-shaped work envelope that resembles an upside-down shield.
        It has two singularities, one at each of the two lateral "tips" of the shield.

        Because the free travel zone cannot be rationalized by a simple rectangle, we limit travel in two ways:
        1) Do the stepper motors go out of bounds?
        2) Do we get near the lateral X limits specified in the deltaxy config section? This helps to avoid motion singularities.
        """

        start_position, end_position = (
            move.start_pos[:3],
            move.end_pos[:3],
        )  # just X, Y, Z
        x_position = end_position[0]

        left_stepper = self.rail_left.get_steppers()[0]
        right_stepper = self.rail_right.get_steppers()[0]
        z_stepper = self.rail_z.get_steppers()[0]

        left_stepper_pos = left_stepper.calc_position_from_coord(end_position)
        right_stepper_pos = right_stepper.calc_position_from_coord(end_position)
        z_stepper_pos = z_stepper.calc_position_from_coord(end_position)

        left_stepper_out_of_bounds = (
            left_stepper_pos > self.left_rail_max_position_mm
        ) or (left_stepper_pos < self.left_rail_min_position_mm)
        right_stepper_out_of_bounds = (
            right_stepper_pos > self.right_rail_max_position_mm
        ) or (right_stepper_pos < self.right_rail_min_position_mm)
        z_stepper_out_of_bounds = (z_stepper_pos > self.z_rail_max_position_mm) or (
            z_stepper_pos < self.z_rail_min_position_mm
        )
        tool_x_out_of_bounds = (
            x_position > self.toolhead_limit_max_x_mm + self.toolhead_offset_x_mm
        ) or (x_position < self.toolhead_limit_min_x_mm + self.toolhead_offset_x_mm)

        if (
            left_stepper_out_of_bounds
            or right_stepper_out_of_bounds
            or z_stepper_out_of_bounds
            or tool_x_out_of_bounds
        ):
            raise move.move_error("Move Out of Bounds!")

        # Limit Z Velocity
        move_distance = move.move_d
        z_distance = abs(move.axes_d[2])
        z_ratio = z_distance / move_distance
        if z_ratio:
            move.limit_speed(self.max_velocity_z / z_ratio, self.max_accel_z / z_ratio)

    def get_status(self, eventtime):

        axes = ""
        if self.homed_axes[0] or self.homed_axes[1]:
            axes += "xy"
        if self.homed_axes[2]:
            axes += "z"

        return {
            "homed_axes": axes,
            # 'axis_minimum': self.axes_min,
            # 'axis_maximum': self.axes_max,
        }


def load_kinematics(toolhead, config):
    return DeltaXYKinematics(toolhead, config)
