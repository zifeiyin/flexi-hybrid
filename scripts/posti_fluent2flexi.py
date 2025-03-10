import numpy as np
import pandas as pd
import h5py
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


def fluent_to_flexi(fluentfile, meshfile, statefile):
    df = pd.read_csv(fluentfile, engine="python", sep=r",\s*")
    kdt = scipy.spatial.KDTree(df[["x-coordinate", "y-coordinate"]])
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
                    distance, index = kdt.query(x[i, j, k, 0:2])
                    max_distance = max(distance, max_distance)
                    rho = df.loc[index, "density"]
                    u = df.loc[index, "x-velocity"]
                    v = df.loc[index, "y-velocity"]
                    w = 0.0
                    p = df.loc[index, "pressure"]
                    tke = max(df.loc[index, "turb-kinetic-energy"], 1.0e-16)
                    omg = max(df.loc[index, "specific-diss-rate"], 1.0e-16)
                    rhoE = (
                        p / (1.4 - 1.0) + 0.5 * rho * (u**2 + v**2 + w**2) + rho * tke
                    )
                    g = np.sqrt(1.0 / (0.09 * omg))
                    U[iElem, i, j, k] = np.array(
                        [rho, rho * u, rho * v, rho * w, rhoE, rho * tke, rho * g]
                    )
        print(f"{iElem}/{nElems}", end="\r")
    print(f"Max distance: {max_distance}")
    state["DG_Solution"][:, :, :, :, :] = U
    mesh.close()
    state.close()


fluent_to_flexi(
    "./rae2822_fluent_refined.csv",
    "./rae2822_mesh.h5",
    "./rae2822_State_0000000.000000000.h5",
)
