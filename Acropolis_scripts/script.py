from acropolis.moreModels import ourDecayModel
import numpy as np

def main():
    res=np.empty((0,29))
    for logEps in np.arange(-4.,-0.5,0.5):
        for logTanT in np.arange(-5.,-0.5,0.5):
            print("Starting caclulation")
            #ourDecayModel(MD1,MD2,MZp,tan_theta,epsilon,alpha_X)
            model = ourDecayModel(500.,525.,1500.,10.**logTanT,10.**logEps,1/(4*np.pi))
            if model._tau() < 1e3:
                continue
            Yf = model.run_disintegration()
            res = np.append(res,np.array([np.concatenate(([logEps, logTanT], Yf.flatten()))]),axis=0)
            print("Ending calculation")
    f=open("MD1_500_MD2_525_MZp_300.dat","w")
    np.savetxt(f,res)
    f.close()


if __name__ == "__main__":
    main()
