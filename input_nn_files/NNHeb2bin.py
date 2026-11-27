"""
This program converts *.dat files as they are present in the /data/nas/ME/NN
sub-directories of strongint at TU-Darmstadt
into the input file format expected by NuHamil-public.

Date: 261123
Author Lyding Anrie Brumm
"""

from functools import reduce
from pathlib import Path
from argparse import ArgumentParser
from typing import BinaryIO, Literal, cast
from dataclasses import dataclass
import numpy as np

# ==============================================================================
# This unit conversion factor converts from the pi/2 * HBARC * M units used
# by Kai Hebeler`s *.dat files to the (MeV)**2 * (fm)**3 units used by
# Takayuki's NuHamil
M = 938.27
HBARC = 197.327
UNIT_CONVERSION_FACTOR = 2 * (HBARC) ** 2 / (np.pi * M)
# ==============================================================================


@dataclass
class Arguments:
    """
    For typed arguments
    """

    Jmax: int
    fn_format: str
    use_Vc: bool
    Np: int
    out: str
    path: Path


def get_kmesh(fn: Path, Np: int) -> list[list[str]]:
    """
    Read the kmesh of a *.dat file
    """
    with open(fn) as fd:
        result = list(map(lambda e: e.split(), fd.read().split("\n")[:-1][:Np]))
    return result


def compare_kmesh(fn: Path, kmesh1: list[list[str]], Np: int) -> bool:
    """
    Compare the kmesh of a *.dat file with an already read kmesh
    """
    kmesh2 = get_kmesh(fn, Np)
    return kmesh2 == kmesh1


def sanity_check_data_files(fns: list[Path], Np: int):
    """
    Perform some sanity checks on all data files within a directory

    For now sanity is declared if:
        1. All *.dat files come with the same kmesh
    """
    kmesh = None
    for fn in fns:
        if kmesh is None:
            kmesh = get_kmesh(Path(fn), Np)
        else:
            if compare_kmesh(Path(fn), kmesh, Np) is False:
                return False, f"{fn} is incompatible k mesh"
    return True, ""


def generate_list_of_channels(
    Jmax: int,
) -> tuple[list[list[dict[str, str]]], list[dict[str, str]]]:
    """
    Generate a list of channels (and their
    """
    singlets: list[list[dict[str, str]]] = []
    triplets: list[list[dict[str, str]]] = []
    for J in range(Jmax + 1):
        for T in [0, 1]:
            for t in range(-T, T + 1):
                # This has to be done prophylactically to group
                # all three states making up the triplet into
                # a sub-list
                # This might lead to empty lists however which
                # are removed at the end of the function
                singlets.append([])
                triplets.append([])
                for S in [0, 1]:
                    for s1 in range(-S, S + 1):
                        for s2 in range(-S, S + 1):
                            singlet = True
                            L1 = J + s1
                            L2 = J + s2
                            if L1 < 0 or L2 < 0:
                                continue
                            if (L1 + S + T) % 2 == 0:
                                continue
                            if (L2 + S + T) % 2 == 0:
                                continue
                            if not (abs(L1 - S) <= J and J <= L1 + S):
                                continue
                            if not (abs(L2 - S) <= J and J <= L2 + S):
                                continue
                            if s1 < 0 and s2 < 0:
                                subtype = "--"
                                singlet = False
                            elif s1 > 0 and s2 > 0:
                                subtype = "++"
                                singlet = False
                            elif s1 < 0 and s2 > 0:
                                subtype = "-+"
                                singlet = False
                            elif s1 > 0 and s2 < 0:
                                subtype = "+-"
                                singlet = False
                            elif s1 == 0 and s2 == 0:
                                subtype = "=="
                                singlet = True
                            else:
                                raise RuntimeError(f"Should never happend: {s1} {s2}")
                            if singlet:
                                singlets[-1].append(
                                    {
                                        "S": str(S),
                                        "L1": str(L1),
                                        "L2": str(L2),
                                        "J": str(J),
                                        "T": str(T),
                                        "Tz": str(t),
                                        "sub": "==",
                                    }
                                )
                            else:
                                triplets[-1].append(
                                    {
                                        "S": str(S),
                                        "L1": str(L1),
                                        "L2": str(L2),
                                        "J": str(J),
                                        "T": str(T),
                                        "Tz": str(t),
                                        "sub": subtype,
                                    }
                                )
    singlets = list(filter(lambda e: len(e) > 0, singlets))
    triplets = list(filter(lambda e: len(e) > 0, triplets))
    result = singlets + triplets
    return result, reduce(lambda a, b: a + b, result)


def channel2fn(
    path: Path, fn_format: str, channel: dict[str, str], use_Vc: bool = False
) -> Path:
    """
    Uses a file name format and a channel tuple to generate a file name
    """
    Tz_str = {-1: "nn", 0: "np", 1: "pp"}[int(channel["Tz"])]
    if use_Vc and channel["Tz"] == "pp":
        return path / Path(
            fn_format.format(
                S=channel["S"],
                L1=channel["L1"],
                L2=channel["L2"],
                J=channel["J"],
                T=channel["T"],
                Tz=Tz_str + "_Vc",
            )
        )
    else:
        return path / Path(
            fn_format.format(
                S=channel["S"],
                L1=channel["L1"],
                L2=channel["L2"],
                J=channel["J"],
                T=channel["T"],
                Tz=Tz_str,
            )
        )


def read_matrix_elements(fn: str, Np: int):
    """
    Read the matrix elements of a data file

    The file format is:

    k_weight k
        ...
    k1 k2 matrix_element
        ...

    """
    with open(fn, "r") as fd:
        # Skip the header and read the third line of the file
        data = list(map(lambda e: float(e.split()[2]), fd.read().split("\n")[Np:-1]))
        return np.array(data, dtype=np.float64)


def write_main_header(
    fd: BinaryIO,
    NMesh: int,
    Jmax: int,
    NChan: int,
    ks: list[float],
    k_weights: list[float],
    byte_order: Literal["little", "big"] = "little",
):
    _ = fd.write(NMesh.to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(Jmax.to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(NChan.to_bytes(4, byteorder=byte_order, signed=True))
    ks2 = np.asarray(ks, dtype=np.float64)
    k_weights2 = np.asarray(k_weights, dtype=np.float64)
    _ = fd.write(ks2.tobytes())
    _ = fd.write(k_weights2.tobytes())


def write_channel(
    fd: BinaryIO,
    path: Path,
    fn_format: str,
    use_Vc: bool,
    Ndim: int,
    channel: list[dict[str, str]],
    byte_order: Literal["little", "big"] = "little",
):
    # Only one sub-channel
    if len(channel) == 1:
        channel_name = channel2fn(path, fn_format, channel[0], use_Vc)
        data = read_matrix_elements(str(channel_name), Ndim)
        data = np.asarray(data, dtype=np.float64)
        full_Ndim = Ndim
    # --, -+, +-, ++ case
    elif len(channel) == 4:
        datas = [np.array([]), np.array([]), np.array([]), np.array([])]
        for i in range(4):
            channel_name = channel2fn(path, fn_format, channel[i], use_Vc)
            tmp_data = read_matrix_elements(str(channel_name), Ndim)
            tmp_data = np.asarray(tmp_data, dtype=np.float64).reshape((Ndim, Ndim))
            datas[{"--": 0, "-+": 1, "+-": 2, "++": 3, "==": 4}[channel[i]["sub"]]] = (
                tmp_data
            )
        full_Ndim = 2 * Ndim
        data = np.block([[datas[0], datas[1]], [datas[2], datas[3]]]).flatten()
        data = np.asarray(data, dtype=np.float64)
    else:
        raise RuntimeError("Should never happen!")

    data = UNIT_CONVERSION_FACTOR * data

    P: int = 1 if int(channel[0]["L1"]) % 2 == 0 else -1
    _ = fd.write(int(channel[0]["J"]).to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(P.to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(int(channel[0]["S"]).to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(int(channel[0]["Tz"]).to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(full_Ndim.to_bytes(4, byteorder=byte_order, signed=True))
    _ = fd.write(data.tobytes())


if __name__ == "__main__":
    parser = ArgumentParser()
    _ = parser.add_argument("--Jmax", type=int, required=True)
    _ = parser.add_argument("--fn_format", type=str, required=True)
    _ = parser.add_argument("--use_Vc", default=False, action="store_true")
    _ = parser.add_argument("--Np", required=True, type=int)
    _ = parser.add_argument("--out", required=True, type=str)
    _ = parser.add_argument("path", type=Path)
    args = parser.parse_args()
    args = Arguments(
        Jmax=cast(int, args.Jmax),
        fn_format=cast(str, args.fn_format),
        use_Vc=cast(bool, args.use_Vc),
        Np=cast(int, args.Np),
        out=cast(str, args.out),
        path=cast(Path, args.path),
    )

    print("Generating channels.")
    channels, channels_flattened = generate_list_of_channels(args.Jmax)
    print("Generating channels file names.")
    channel_names = list(
        map(
            lambda channel: channel2fn(args.path, args.fn_format, channel, args.use_Vc),
            channels_flattened,
        )
    )
    print("Sanity checking channel files.")
    if sanity_check_data_files(channel_names, args.Np):
        print("Channel files are found to be sane.")
    else:
        print("Channel files are not found to be sane.")
        exit(-1)
    print(f"Reading mesh from {str(channel_names[0])}")
    mesh_str = get_kmesh(channel_names[0], args.Np)
    ks = list(map(lambda e: float(e[1]), mesh_str))
    k_weights = list(map(lambda e: float(e[0]), mesh_str))
    print(f"Opening {args.out} for (binary) writing.")
    with open(args.out, "bw") as fd:
        print("Writing main header")
        write_main_header(fd, args.Np, args.Jmax, len(channels), ks, k_weights)
        for channel in channels:
            write_channel(
                fd,
                args.path,
                args.fn_format,
                args.use_Vc,
                args.Np,
                channel,
                byte_order="little",
            )
    print(f"Closing {args.out} for (binary) writing.")
