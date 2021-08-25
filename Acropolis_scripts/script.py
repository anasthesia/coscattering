from acropolis.moreModels import ourDecayModel
import numpy as np

def main():
    res=np.empty((0,29))
    for logEps in np.arange(-4.,-0.5,0.5):
        for logTanT in np.arange(-5.,-0.5,0.5):
            model = ourDecayModel(300.,315.,900.,10.**logTanT,10.**logEps,1/(4*np.pi))
            if model._tau() < 1e3:
                continue
            Yf = model.run_disintegration()
            res = np.append(res,np.array([np.concatenate(([logEps, logTanT], Yf.flatten()))]),axis=0)

    f=open("MD1_300_MD2_315_MZp_900.dat","w")
    np.savetxt(f,res)
    f.close()


if __name__ == "__main__":
    main()
