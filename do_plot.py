from matplotlib import pyplot as pl
import seaborn as sb

def denser(data,label=None,color=None,x1=None,x2=None,shade=None,save=None):
    
    # Important parameters
    # color
    # xlim
    # shade
    # opacity
    
    
    sb.kdeplot(data,label=label,color=color,shade=shade)
    pl.legend()
    pl.xlim(x1,x2)
    
    pl.yticks([])



    pl.xlabel('')

    pl.ylabel('')

    
    if save:
        
        pl.savefig('denser_plot',dpi=300)