import os
os.environ["PATH"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'
from diagrams import Diagram
from diagrams import Cluster
from diagrams.aws.compute import *
from diagrams.aws.database import *
from diagrams.aws.network import ELB, Route53
from diagrams.aws.analytics import *
from diagrams.aws.business import *

with Diagram("Particle Collision Test", show=False): # Keep
    collision = ComputeOptimizer("collision") # Keep
    vphoton = Workmail("virtual photon") # Keep
    with Cluster("Initial Particles"): # Keep
        proton = A4B("proton")  # Keep
        electron = Chime("electron") # Keep
    with Cluster("Proton"): # Keep
        quark = AMI("Down")

        double_quark = AMI("Up, Up")
    Lund = Compute("Lund")
    with Cluster("Post Collision"):
        struck_quark = AMI("Struck Down")

        diquark = Compute("diquark 2101")

        mhadron0 = Aurora("Rho+")

    with Cluster("Final State"):
        with Cluster("Hadrons"):
            hadron0 = Compute("Pi+")

            hadron1 = Compute("Pi0")

            hadron2 = Compute("Pi-")

            hadron3 = Compute("Proton")

            hadron4 = Compute("Rho0")

            hadron5 = Compute("Pi-")

            hadron6 = Compute("Pi+")

        with Cluster("Decay Photon Pairs"):
            fphoton0 = Neptune("Photon")

            fphoton1 = Neptune("Photon")

        scattered_electron = Chime("Scattered Electron")
    hadron1 >> fphoton1
    fphoton1 - fphoton0
    hadron4 >> hadron5
    hadron4 >> hadron6
    proton >> double_quark
    proton >> quark >> collision >> struck_quark >> Lund 
    double_quark >> mhadron0 >> Lund
    double_quark >> diquark >> Lund
    Lund >> hadron2
    Lund >> hadron3
    Lund >> hadron4
    electron >> vphoton >> collision
    electron >> scattered_electron
