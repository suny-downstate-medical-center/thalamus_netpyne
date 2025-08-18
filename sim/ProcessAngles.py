class AngleProcessor:
    def __init__(self):
        self.reference_angles_ = list(range(0, 360, 45))
        self.store_reference_angles_ = list(range(0, 360, 45))
        self.reference_angles = self._convert_odd_to_even(self.reference_angles_)
        self.angle_slices_ = self._unpack_even_angles(self.reference_angles)
        self.angle_slices = self._fix_negative_values(self.angle_slices_)
        self.angle_slices_relative_y, self.angle_slices_relative_y_inverted = self._create_relative_value_dicts()

    def _convert_odd_to_even(self, angles):
        reference_angles = []
        for angle in angles:
            if angle % 2 != 0:
                angle_ = angle - 1
            else:
                angle_ = angle
            reference_angles.append(angle_)
        return reference_angles

    def _unpack_even_angles(self, angles):
        angle_slices_ = {}
        for angle in angles:
            angle_slices_.update({angle: [i for i in range(angle - 22, angle + 23, 2)]})
        return angle_slices_

    def _fix_negative_values(self, angle_slices_):
        angle_slices = {}
        for angle in angle_slices_.keys():
            angle_slices.update({angle: []})
            for angle_val_ in angle_slices_[angle]:
                if angle_val_ < 0:
                    angle_val = angle_val_ + 360
                else:
                    angle_val = angle_val_
                angle_slices[angle].append(angle_val)
        return angle_slices

    def _create_relative_value_dicts(self):
        angle_slices_relative_y = {}
        angle_slices_relative_y_inverted = []
        for angle_ind, angle in enumerate(self.angle_slices.keys()):
            angle_slices_relative_y.update({self.store_reference_angles_[angle_ind]: [val / 360 for val in self.angle_slices[angle]]})
            for val in self.angle_slices[angle]:
                angle_slices_relative_y_inverted.append((val / 360, self.store_reference_angles_[angle_ind]))
        return angle_slices_relative_y, angle_slices_relative_y_inverted

    def get_closest_angle(self, value):
        if not 0 <= value <= 1:
            raise ValueError("Input value must be between 0 and 1")
        closest_value = min(self.angle_slices_relative_y_inverted, key=lambda x: abs(x[0] - value))
        return closest_value[1]

# Usage
angle_processor = AngleProcessor()
print(angle_processor.reference_angles)
print(angle_processor.angle_slices)
print(angle_processor.angle_slices_relative_y)
print(angle_processor.angle_slices_relative_y_inverted)
print(angle_processor.get_closest_angle(0.25))
