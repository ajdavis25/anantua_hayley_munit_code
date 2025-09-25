import os, tempfile, pytest, numpy as np
from pathlib import Path
from src.plot.creating_images import infer_positron_ratio, plotPositronTestFrame


def test_infer_positron_ratio_from_filename():
    assert infer_positron_ratio("model_pos0.h5") == 0.0
    assert infer_positron_ratio("model_pos1.h5") == 1.0
    assert infer_positron_ratio("model_pos0.5.h5") == 0.5
    assert infer_positron_ratio("model_nopos.h5") is None


@pytest.fixture
def dummy_h5_file(tmp_path):
    # create a minimal dummy .h5 file with necessary structure
    import h5py
    path = tmp_path / "dummy_0.h5"
    with h5py.File(path, "w") as f:
        header = f.create_group("header")
        camera = header.create_group("camera")
        camera.create_dataset("dx", data=1.0)
        camera.create_dataset("dy", data=1.0)
        header.create_dataset("dsource", data=1.0)
        units = header.create_group("units")
        units.create_dataset("L_unit", data=1.0)
        header.create_dataset("scale", data=1.0)
        f.create_dataset("unpol", data=np.zeros((2, 2)))
        f.create_dataset("pol", data=np.zeros((2, 2, 4)))
    return path


def test_plotPositronTestFrame_runs_without_error(dummy_h5_file):
    output_png = dummy_h5_file.with_suffix(".png")
    plotPositronTestFrame(str(dummy_h5_file), output=str(output_png))
    assert output_png.exists()
