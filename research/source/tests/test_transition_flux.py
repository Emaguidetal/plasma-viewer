import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from rf_plasma.transition_rf import _field_flux, _power_wave_ratio
from rf_plasma.validate_results import _reciprocity_error


class TransitionFluxTests(unittest.TestCase):
    def test_power_wave_uses_each_mode_impedance_and_rejects_invalid_inputs(self):
        incident = np.array([2.0 + 0j])
        outgoing = np.array([1.0j])
        result = _power_wave_ratio(outgoing, incident, np.array([800.0]), np.array([200.0]))
        self.assertAlmostEqual(abs(result[0]) ** 2, 1.0)
        with self.assertRaises(ValueError):
            _power_wave_ratio(outgoing, np.array([0j]), np.array([800.0]), np.array([200.0]))
        with self.assertRaises(ValueError):
            _power_wave_ratio(outgoing, incident, np.array([800.0]), np.array([200.0 + 1j]))

    def test_field_flux_integrates_signed_poynting_vector(self):
        x = np.array([0.0, 0.5, 2.0])
        y = np.array([0.0, 1.3, 3.0])
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            for name in ("electric", "magnetic"):
                with h5py.File(output / f"input_{name}_flux.h5", "w") as handle:
                    mesh = handle.create_group("Mesh")
                    mesh.create_dataset("x", data=x)
                    mesh.create_dataset("y", data=y)
                    mesh.create_dataset("z", data=np.array([0.0]))
                    data = np.zeros((3, len(x), len(y), 1), dtype=complex)
                    if name == "electric":
                        data[1] = 1j
                    else:
                        data[0] = -0.25j
                    dataset = handle.create_dataset("FieldData/FD/f0", data=data)
                    dataset.attrs["frequency"] = 3e9
            self.assertAlmostEqual(_field_flux(output, "input", 3e9), 0.75)
            with h5py.File(output / "input_magnetic_flux.h5", "r+") as handle:
                handle["FieldData/FD/f0"][0] = 0.25j
            self.assertAlmostEqual(_field_flux(output, "input", 3e9), -0.75)
            with self.assertRaises(RuntimeError):
                _field_flux(output, "input", 2e9)
            with self.assertRaises(ValueError):
                _field_flux(output, "input", float("nan"))
            with h5py.File(output / "input_magnetic_flux.h5", "r+") as handle:
                handle["Mesh/y"][0] = -1
            with self.assertRaises(RuntimeError):
                _field_flux(output, "input", 3e9)

    def test_reciprocity_checks_complex_phase_and_frequency(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            forward, reverse = root / "forward", root / "reverse"
            forward.mkdir()
            reverse.mkdir()
            np.savez(forward / "modal-data.npz", frequency_hz=[3e9], **{"RECT-TE10_raw_power_wave": [1 + 0j]})
            np.savez(reverse / "modal-data.npz", frequency_hz=[3e9], **{"TE11-90_raw_power_wave": [1j]})
            self.assertAlmostEqual(_reciprocity_error(forward, reverse, 3e9), np.sqrt(2))
            np.savez(reverse / "modal-data.npz", frequency_hz=[2e9], **{"TE11-90_raw_power_wave": [1j]})
            with self.assertRaises(ValueError):
                _reciprocity_error(forward, reverse, 3e9)


if __name__ == "__main__":
    unittest.main()
