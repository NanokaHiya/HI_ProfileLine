#!/usr/bin/env python
#imported. minor changes.
from optparse import OptionParser
import numpy as np
import numpy.linalg as LA


def get_shape(x, mpart, Rb=20., decrease=True):
    # Rb=20kpc, within which the shape is calcualted
    s = 1.
    q = 1.

    Tiv = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])
    # tmp = Tiv

    order = [2, 1, 0]

    Vei = np.zeros((3, 3))

    dq = 10000.
    ds = 10000.

    while (dq > 0.01 or ds > 0.01):
        # in eigenvector coordinates
        y = np.transpose(np.dot(Tiv, np.transpose(x)))
        rn0 = np.sqrt(np.power(y[:, order[2]], 2.) +
                      np.power(y[:, order[1]], 2.)/q/q +
                      np.power(y[:, order[0]], 2.)/s/s)
        print('rn0={}'.format(rn0))
        ind = np.where(rn0 < Rb)[0]
        print('ind = ')
        print(ind)
        # Np = ind.shape[0]

        y1 = y[ind, 0]
        print(y1)
        y2 = y[ind, 1]
        y3 = y[ind, 2]
        rn = rn0[ind]
        mass=mpart[ind]
        mtot=np.sum(mass)
        print('mtot: {}'.format(mtot))         
        I11 = np.sum(y1*y1/np.power(rn, 2)*mass)
        I22 = np.sum(y2*y2/np.power(rn, 2)*mass)
        I33 = np.sum(y3*y3/np.power(rn, 2)*mass)
        I12 = np.sum(y1*y2/np.power(rn, 2)*mass)
        I13 = np.sum(y1*y3/np.power(rn, 2)*mass)
        I23 = np.sum(y2*y3/np.power(rn, 2)*mass)

        II = [[I11, I12, I13],
              [I12, I22, I23],
              [I13, I23, I33]]

        print(II)
        II=II/mtot

        D, A = LA.eig(II)
        # print 'eigenvalues'
        # print D
        # print  'eigenvectors'
        # print A
        order = np.argsort(D)  # a=order2,b=order1,c=order0
        la = np.sqrt(D[order[2]])
        lb = np.sqrt(D[order[1]])
        lc = np.sqrt(D[order[0]])

        dq = np.abs(q-lb/la)
        ds = np.abs(s-lc/la)

        q = lb/la
        s = lc/la

        Tiv = np.dot(LA.inv(A), Tiv)

    # rba = q
    # rca = s
    if decrease:
        Tiv = Tiv[order[::-1], :]
    else:
        Tiv = Tiv[order, :]
    # right hand coordinate
    if np.dot(Tiv[2, :], np.cross(Tiv[0, :], Tiv[1, :])) < 0:
        Tiv[1, :] *= -1
    # eigen vectors Vei[:,0] Vei[:,1] Vei[:,2]
    Vei = LA.inv(Tiv)

    d = np.array([0, 0, 1])
    costh = np.dot(Vei[:, 2], d) / np.sqrt(np.dot(Vei[:, 2], Vei[:, 2])) /\
        np.sqrt(np.dot(d, d))
    # angle between longest axis (z' direction) and LOS (i.e. z direction)
    angle = np.arccos(costh)*180./np.pi

    ba = q
    ca = s
    # print ' b/a= {:.2f}'.format(q),'  c/a = {:.2f}'.format(s)
    # print "rotation angle=",angle

    return ba, ca, angle, Tiv

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-r', action='store', type='float', dest='R',
                      default=10.0, help='calculation radius')
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print('Error - please provide a folder name')
        exit(1)
    data = np.load('{}/part_org.npy'.format(args[0]))  # data without rotation
    xstar = data[0]  # star particle position, N*3 array [kpc]
    vstar = data[1]  # star particle velocity, N*3 array [km/s]
    mstar = data[2]  # star particle mass, N*1 array [10^10 M_solar]
    xdark = data[3]  # dark matter particle position, N*3 array [kpc]
    vdark = data[4]  # dark matter particle velocity, N*3 array [km/s]
    mdark = data[5]  # dark matter particle mass, N*1 array [10^10 M_solar]
    print('Read data success, please wait for a few minutes for axis '
          'ratio calculation')
    ba, ca, angle, Tiv = get_shape(xstar, mstar, Rb=options.R)
    xstar_axis = np.dot(Tiv, xstar.T).T
    vstar_axis = np.dot(Tiv, vstar.T).T

    xdark_axis = np.dot(Tiv, xdark.T).T
    vdark_axis = np.dot(Tiv, vdark.T).T

    print('Calculation radius: {:.2f}'.format(options.R))
    print('Axis ratios b/a={:.2f}  c/a={:.2f}'.format(ba, ca))
    print('Rotation matrix:')
    print('{:6.3f} {:6.3f} {:6.3f}'.format(Tiv[0, 0], Tiv[0, 1], Tiv[0, 2]))
    print('{:6.3f} {:6.3f} {:6.3f}'.format(Tiv[1, 0], Tiv[1, 1], Tiv[1, 2]))
    print('{:6.3f} {:6.3f} {:6.3f}'.format(Tiv[2, 0], Tiv[2, 1], Tiv[2, 2]))

