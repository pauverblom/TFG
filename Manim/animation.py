from manim import *

import numpy as np

class FixedParticleStep(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75*DEGREES, theta=-45*DEGREES)
        max_angle = PI/3
        num_points = 100

        # Initialize axes with manual labels
        axes = ThreeDAxes(
            x_range=[-1, 1, 1],
            y_range=[-1, 1, 1],
            z_range=[-1, 1, 1],
            axis_config={"color": GREY}
        )
        # Add axis labels manually
        x_label = axes.get_x_axis_label("x", edge=RIGHT, direction=RIGHT)
        y_label = axes.get_y_axis_label("y", edge=UP, direction=UP)
        z_label = axes.get_z_axis_label("z", edge=OUT, direction=OUT)
        self.add(axes, x_label, y_label, z_label)

        # Create initial random direction
        last_dir = normalize(np.random.randn(3))
        last_arrow = Arrow3D(
            start=ORIGIN,
            end=last_dir,
            color=RED,
            resolution=8
        )
        
        # Label positioning
        last_label = Text("Last direction", font_size=16, color = RED)
        self.add_fixed_in_frame_mobjects(last_label)
        last_label.next_to(last_arrow.get_end(), UP+LEFT, buff=0)

        # Create spherical cap cone
        cone_height = np.cos(max_angle)
        cone_radius = np.sin(max_angle)
        cap_cone = Cone(
            direction=IN,
            base_radius=cone_radius,
            height=cone_height,
            fill_color=WHITE, 
            fill_opacity=0.1,  # Low opacity
            stroke_opacity=0
        )


        rotation_axis = normalize(np.cross(OUT, last_dir))
        rotation_angle = np.arccos(np.dot(OUT, last_dir))
        
        # Generate random points in spherical cap
        points = VGroup()
        for _ in range(num_points):
            theta = np.arccos(1 - (1 - np.cos(max_angle)) * np.random.rand())
            phi = 2 * np.pi * np.random.rand()
            point = Dot3D(
            point=spherical_to_cartesian(theta, phi, 1),
            color=BLUE_B,
            radius=0.02
            )
            points.add(point)

        # Add initial elements
        self.add(last_arrow, last_label, cap_cone, points)
        self.wait(1)

        # Rest of the animation remains the same...
        # Select and highlight random point
        selected_point = points[np.random.randint(len(points))].copy()
        selected_point.set_color(YELLOW)
        self.play(
            FadeOut(points),
            selected_point.animate.scale(2),
            run_time=1
        )
        self.wait(0.5)

        # Rotate selected point and cone to last direction
        self.play(
            Rotate(
                selected_point,
                angle=rotation_angle,
                axis=rotation_axis,
                about_point=ORIGIN
            ),
            Rotate(
                cap_cone,
                angle=rotation_angle,
                axis=rotation_axis,
                about_point=ORIGIN
            ),
            run_time=2
        )
        self.wait(1)

        # Create new direction vector and label
        new_dir = normalize(selected_point.get_center())
        new_arrow = Arrow3D(
            start=ORIGIN,
            end=new_dir,
            color=YELLOW,
            resolution=8
        )
        new_label = Text("New direction", font_size=16, color=YELLOW)
        self.add_fixed_in_frame_mobjects(new_label)
        new_label.next_to(new_arrow.get_end(), DOWN+RIGHT, buff=0)

        # Animate transition
        self.play(
            ReplacementTransform(last_arrow, new_arrow),
            FadeIn(new_label),
            FadeOut(selected_point, cap_cone, last_label),
            run_time=1.5
        )
        self.wait(2)

def spherical_to_cartesian(theta, phi, r=1):
    return np.array([
        r * np.sin(theta) * np.cos(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(theta)
    ])
