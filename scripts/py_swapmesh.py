#!/usr/bin/env python3

### created by Sijie Wang on March 6, 2025 ###

import argparse

import h5py
import numpy as np
import scipy.spatial
import scipy.special


def vandermonde(input, output):
    m = len(output)
    n = len(input)
    A = np.ones((m, n))
    for i in range(m):
        for j in range(n):
            for l in range(n):
                if l != j:
                    A[i, j] *= (output[i] - input[l]) / (input[j] - input[l])
    return A


def Ngeo_to_N_vdm(Ngeo, N):
    x0 = np.array([-1 + i * 2 / Ngeo for i in range(Ngeo + 1)])
    x1 = np.array([np.cos(i / Ngeo * np.pi) for i in range(Ngeo + 1)])
    x2 = np.array([np.cos(i / N * np.pi) for i in range(N + 1)])
    x3, _ = scipy.special.roots_legendre(N + 1)
    A1 = vandermonde(x0, x1)
    A2 = vandermonde(x1, x2)
    A3 = vandermonde(x2, x3)
    return A3 @ A2 @ A1


def change_basis(U, A):
    N1, N2, N3, Nv = U.shape
    No, Ni = A.shape
    assert N1 == N2 == N3 == Ni
    U1 = np.zeros((Ni, Ni, No, Nv))
    U2 = np.zeros((Ni, No, No, Nv))
    U3 = np.zeros((No, No, No, Nv))
    for i in range(Ni):
        for j in range(Ni):
            for k in range(No):
                for l in range(Ni):
                    U1[i, j, k] += A[k, l] * U[i, j, l]
    for i in range(Ni):
        for j in range(No):
            for l in range(Ni):
                for k in range(No):
                    U2[i, j, k] += A[j, l] * U1[i, l, k]
    for i in range(No):
        for l in range(Ni):
            for j in range(No):
                for k in range(No):
                    U3[i, j, k] += A[i, l] * U2[l, j, k]
    return U3


def buildCoords(NodeCoords, nElems, Ngeo, N):
    x = np.zeros((nElems, N + 1, N + 1, N + 1, 3))
    A = Ngeo_to_N_vdm(Ngeo, N)
    iNode = 0
    for iElem in range(nElems):
        x_geo = np.zeros((Ngeo + 1, Ngeo + 1, Ngeo + 1, 3))
        for i in range(Ngeo + 1):
            for j in range(Ngeo + 1):
                for k in range(Ngeo + 1):
                    x_geo[i, j, k] = NodeCoords[iNode]
                    iNode += 1
        x[iElem] = change_basis(x_geo, A)
    return x


def construct_dataset(meshfile, statefile):
    mesh = h5py.File(meshfile)
    Ngeo = int(mesh.attrs["Ngeo"][0])
    nElems = int(mesh.attrs["nElems"][0])
    NodeCoords = mesh["NodeCoords"]
    state = h5py.File(statefile)
    N = int(state.attrs["N"][0])
    nDoFs_elem = (N + 1) ** 3
    nDoFs = nElems * nDoFs_elem
    x = buildCoords(NodeCoords, nElems, Ngeo, N)
    x.resize((nDoFs, 3))
    U = np.zeros(state["DG_Solution"].shape)
    U[:] = state["DG_Solution"]
    U.resize((nDoFs, 7))
    mesh.close()
    state.close()
    return x, U


def swapmesh2(meshfile_old, statefile_old, meshfile, statefile):
    x_old, U_old = construct_dataset(meshfile_old, statefile_old)
    kdt = scipy.spatial.KDTree(x_old)
    mesh = h5py.File(meshfile)
    Ngeo = int(mesh.attrs["Ngeo"][0])
    nElems = int(mesh.attrs["nElems"][0])
    NodeCoords = mesh["NodeCoords"]
    state = h5py.File(statefile, "r+")
    N = int(state.attrs["N"][0])
    U = np.zeros(state["DG_Solution"].shape)
    A = Ngeo_to_N_vdm(Ngeo, N)
    iNode = 0
    max_distance = 0.0
    for iElem in range(nElems):
        x_geo = np.zeros((Ngeo + 1, Ngeo + 1, Ngeo + 1, 3))
        for i in range(Ngeo + 1):
            for j in range(Ngeo + 1):
                for k in range(Ngeo + 1):
                    x_geo[i, j, k] = NodeCoords[iNode]
                    iNode += 1
        x = change_basis(x_geo, A)
        for i in range(N + 1):
            for j in range(N + 1):
                for k in range(N + 1):
                    distance, index = kdt.query(x[i, j, k])
                    max_distance = max(distance, max_distance)
                    U[iElem, i, j, k] = U_old[index]
        print(f"{iElem}/{nElems}", end="\r")
    print(f"Max distance: {max_distance}")
    state["DG_Solution"][:, :, :, :, :] = U
    mesh.close()
    state.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("meshfile_old")
    parser.add_argument("statefile_old")
    parser.add_argument("meshfile")
    parser.add_argument("statefile")
    args = parser.parse_args()
    swapmesh2(args.meshfile_old, args.statefile_old, args.meshfile, args.statefile)
