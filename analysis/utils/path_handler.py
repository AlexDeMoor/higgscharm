import os
import pathlib


class Paths:

    def __init__(self, eos: bool) -> None:
        if eos:
            # finds the /eos user directory
            user = os.environ["USER"]
            self.root_path = pathlib.Path(f"/eos/user/{user[0]}/{user}/higgscharm")
        else:
            # repository root (<repo>/outputs/ is where run_postprocess.py and the
            # scripts read from; the old value pointed to <repo>/analysis/outputs)
            self.root_path = pathlib.Path(__file__).resolve().parents[2]

    def workflow_path(
        self,
        workflow: str,
        year: str,
        dataset: str
    ) -> pathlib.Path:
        """
        Safely return a path by creating the parent directories to avoid errors when writing to the path.

        Parameters:
        -----------
            path: Path to optionally create and return.
            mkdir: If True, creates the parent directories. If False, it has no effect.

        Returns:
            Input path.
        """
        workflow_path = "/".join(
            [
                elem
                for elem in [
                    workflow,
                    year,
                    dataset
                ]
                if elem is not None
            ]
        )
        # make output directory
        output_path = self.root_path / "outputs" / workflow_path
        if not output_path.exists():
            output_path.mkdir(parents=True)
        return output_path