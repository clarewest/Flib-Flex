#!/usr/bin/python

import Bio
from Bio.PDB import Entity, PDBParser, PPBuilder, Polypeptide
from math import sqrt
from numpy import *
from sys import stdin, argv
import Bio.PDB

def get_angle_list(polypep):
#            """Return the list of dihedral and bond angles."""
        ppl = []
        lng = len(polypep)
        for i in range(0, lng):
                res = polypep[i]
                try:
                        n = res['N'].get_vector()
                        ca = res['CA'].get_vector()
                        c = res['C'].get_vector()
                except Exception:
                        # Some atoms are missing
                        # Phi/Psi cannot be calculated for this residue
                        ppl.append((None, None, None, None, None, None, None, None, None))
                        continue
                # Phi
                if i > 0:
                        rp = polypep[i - 1]
                        try:
                                cp = rp['C'].get_vector()
                                phi = Polypeptide.calc_dihedral(cp, n, ca, c)
                        except Exception:
                                phi = None
                else:
                # No phi for residue 0!
                        phi = None
                # Psi
                if i < (lng - 1):
                        rn = polypep[i + 1]
                        try:
                                nn = rn['N'].get_vector()
                                psi = Polypeptide.calc_dihedral(n, ca, c, nn)
                        except Exception:
                                psi = None
                else:
                # No psi for last residue!
                        psi = None

                # Omega
                if i < (lng-1):
                        rn = polypep[i + 1]
                        try:
                                can = rn['CA'].get_vector()
                                omega = Polypeptide.calc_dihedral(ca, c, nn, can)
                        except Exception:
                                omega = None
                else:
                        omega = None

        #       N-CA-C -- N-CA-C
                # CA-C-N Angle
                if i < (lng -1):
                        cacn = Polypeptide.calc_angle(ca,c,nn)
                else:
                        cacn = None
                # C-N-CA
                if i > 0:
                        cnca = Polypeptide.calc_angle(cp,n,ca)
                else:
                        cnca = None
                # N-CA-C
                ncac = Polypeptide.calc_angle(n,ca,c)
                # CA-C-O
                try:
                        o = res['O'].get_vector()
                except:
                        caco = None
                else:
                        caco = Polypeptide.calc_angle(ca,c,o)
                # Tau
                if i > 1 and i < (lng-1):
                        try:
                                rpp = polypep[i-2]
                                cap = rp['CA'].get_vector()
                                capp = rpp['CA'].get_vector()
                                can = rn['CA'].get_vector()
                        except:
                                tau = None
                        else:
                                tau = Polypeptide.calc_dihedral(capp,cap,ca,can)
                else:
                        tau = None
                # theta
                if i>1 and i < (lng-1):
                        try:
                                can = rn['CA'].get_vector()
                        except:
                                theta = None
                        else:
                                theta = Polypeptide.calc_angle(cap,ca,can)
                else:
                        theta = None
                ppl.append((phi,psi,omega,cacn,cnca,ncac,caco,theta,tau))
        return ppl

if __name__ == '__main__':
        pdb = argv[1]
        cha = argv[2]
        output = open(argv[1]+".angles","w")
        for model in Bio.PDB.PDBParser().get_structure(pdb, pdb+".pdb"):
                chain = model[cha]
                output.write("##### PDB "+pdb+" Chain "+chain.get_id()+"\n" )
                polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides) :
                        j=0
                        chain_angle = []
                        for angles_list in get_angle_list(poly):
                                torsion_angles = []
                                for i in range(len(angles_list)):
                                        try:
                                                torsion_angles.append(float(angles_list[i])*57.2957795)
                                        except:
                                                torsion_angles.append("-999")
                                output.write(str(j+1)+"\t"+Polypeptide.three_to_one(poly[j].resname))
                                
                                for angle in torsion_angles[0:2]:
                                    if angle == "-999":
                                        output.write("\t-999")
                                    else:
                                        output.write("\t{0:.2f}".format(angle))
                                chain_angle.append(torsion_angles)
                                output.write("\n")
                                j=j+1
output.close()
