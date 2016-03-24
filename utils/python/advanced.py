from matplotlib import _cntr as cntr
import matplotlib.gridspec as gridspec
from subprocess import call
from scipy.optimize import curve_fit
from scipy.integrate import cumtrapz
import copy
def load_viridis():
    from matplotlib.colors import LinearSegmentedColormap
    parameters = {'xp': [22.674387857633945, 11.221508276482126, -14.356589454756971, -47.18817758739222, -34.59001004812521, -6.0516291196352654],
                  'yp': [-20.102530541012214, -33.08246073298429, -42.24476439790574, -5.595549738219887, 42.5065445026178, 40.13395157135497],
                  'min_JK': 18.8671875,
                  'max_JK': 92.5}

    cm_data = [[ 0.26700401,  0.00487433,  0.32941519],
           [ 0.26851048,  0.00960483,  0.33542652],
           [ 0.26994384,  0.01462494,  0.34137895],
           [ 0.27130489,  0.01994186,  0.34726862],
           [ 0.27259384,  0.02556309,  0.35309303],
           [ 0.27380934,  0.03149748,  0.35885256],
           [ 0.27495242,  0.03775181,  0.36454323],
           [ 0.27602238,  0.04416723,  0.37016418],
           [ 0.2770184 ,  0.05034437,  0.37571452],
           [ 0.27794143,  0.05632444,  0.38119074],
           [ 0.27879067,  0.06214536,  0.38659204],
           [ 0.2795655 ,  0.06783587,  0.39191723],
           [ 0.28026658,  0.07341724,  0.39716349],
           [ 0.28089358,  0.07890703,  0.40232944],
           [ 0.28144581,  0.0843197 ,  0.40741404],
           [ 0.28192358,  0.08966622,  0.41241521],
           [ 0.28232739,  0.09495545,  0.41733086],
           [ 0.28265633,  0.10019576,  0.42216032],
           [ 0.28291049,  0.10539345,  0.42690202],
           [ 0.28309095,  0.11055307,  0.43155375],
           [ 0.28319704,  0.11567966,  0.43611482],
           [ 0.28322882,  0.12077701,  0.44058404],
           [ 0.28318684,  0.12584799,  0.44496   ],
           [ 0.283072  ,  0.13089477,  0.44924127],
           [ 0.28288389,  0.13592005,  0.45342734],
           [ 0.28262297,  0.14092556,  0.45751726],
           [ 0.28229037,  0.14591233,  0.46150995],
           [ 0.28188676,  0.15088147,  0.46540474],
           [ 0.28141228,  0.15583425,  0.46920128],
           [ 0.28086773,  0.16077132,  0.47289909],
           [ 0.28025468,  0.16569272,  0.47649762],
           [ 0.27957399,  0.17059884,  0.47999675],
           [ 0.27882618,  0.1754902 ,  0.48339654],
           [ 0.27801236,  0.18036684,  0.48669702],
           [ 0.27713437,  0.18522836,  0.48989831],
           [ 0.27619376,  0.19007447,  0.49300074],
           [ 0.27519116,  0.1949054 ,  0.49600488],
           [ 0.27412802,  0.19972086,  0.49891131],
           [ 0.27300596,  0.20452049,  0.50172076],
           [ 0.27182812,  0.20930306,  0.50443413],
           [ 0.27059473,  0.21406899,  0.50705243],
           [ 0.26930756,  0.21881782,  0.50957678],
           [ 0.26796846,  0.22354911,  0.5120084 ],
           [ 0.26657984,  0.2282621 ,  0.5143487 ],
           [ 0.2651445 ,  0.23295593,  0.5165993 ],
           [ 0.2636632 ,  0.23763078,  0.51876163],
           [ 0.26213801,  0.24228619,  0.52083736],
           [ 0.26057103,  0.2469217 ,  0.52282822],
           [ 0.25896451,  0.25153685,  0.52473609],
           [ 0.25732244,  0.2561304 ,  0.52656332],
           [ 0.25564519,  0.26070284,  0.52831152],
           [ 0.25393498,  0.26525384,  0.52998273],
           [ 0.25219404,  0.26978306,  0.53157905],
           [ 0.25042462,  0.27429024,  0.53310261],
           [ 0.24862899,  0.27877509,  0.53455561],
           [ 0.2468114 ,  0.28323662,  0.53594093],
           [ 0.24497208,  0.28767547,  0.53726018],
           [ 0.24311324,  0.29209154,  0.53851561],
           [ 0.24123708,  0.29648471,  0.53970946],
           [ 0.23934575,  0.30085494,  0.54084398],
           [ 0.23744138,  0.30520222,  0.5419214 ],
           [ 0.23552606,  0.30952657,  0.54294396],
           [ 0.23360277,  0.31382773,  0.54391424],
           [ 0.2316735 ,  0.3181058 ,  0.54483444],
           [ 0.22973926,  0.32236127,  0.54570633],
           [ 0.22780192,  0.32659432,  0.546532  ],
           [ 0.2258633 ,  0.33080515,  0.54731353],
           [ 0.22392515,  0.334994  ,  0.54805291],
           [ 0.22198915,  0.33916114,  0.54875211],
           [ 0.22005691,  0.34330688,  0.54941304],
           [ 0.21812995,  0.34743154,  0.55003755],
           [ 0.21620971,  0.35153548,  0.55062743],
           [ 0.21429757,  0.35561907,  0.5511844 ],
           [ 0.21239477,  0.35968273,  0.55171011],
           [ 0.2105031 ,  0.36372671,  0.55220646],
           [ 0.20862342,  0.36775151,  0.55267486],
           [ 0.20675628,  0.37175775,  0.55311653],
           [ 0.20490257,  0.37574589,  0.55353282],
           [ 0.20306309,  0.37971644,  0.55392505],
           [ 0.20123854,  0.38366989,  0.55429441],
           [ 0.1994295 ,  0.38760678,  0.55464205],
           [ 0.1976365 ,  0.39152762,  0.55496905],
           [ 0.19585993,  0.39543297,  0.55527637],
           [ 0.19410009,  0.39932336,  0.55556494],
           [ 0.19235719,  0.40319934,  0.55583559],
           [ 0.19063135,  0.40706148,  0.55608907],
           [ 0.18892259,  0.41091033,  0.55632606],
           [ 0.18723083,  0.41474645,  0.55654717],
           [ 0.18555593,  0.4185704 ,  0.55675292],
           [ 0.18389763,  0.42238275,  0.55694377],
           [ 0.18225561,  0.42618405,  0.5571201 ],
           [ 0.18062949,  0.42997486,  0.55728221],
           [ 0.17901879,  0.43375572,  0.55743035],
           [ 0.17742298,  0.4375272 ,  0.55756466],
           [ 0.17584148,  0.44128981,  0.55768526],
           [ 0.17427363,  0.4450441 ,  0.55779216],
           [ 0.17271876,  0.4487906 ,  0.55788532],
           [ 0.17117615,  0.4525298 ,  0.55796464],
           [ 0.16964573,  0.45626209,  0.55803034],
           [ 0.16812641,  0.45998802,  0.55808199],
           [ 0.1666171 ,  0.46370813,  0.55811913],
           [ 0.16511703,  0.4674229 ,  0.55814141],
           [ 0.16362543,  0.47113278,  0.55814842],
           [ 0.16214155,  0.47483821,  0.55813967],
           [ 0.16066467,  0.47853961,  0.55811466],
           [ 0.15919413,  0.4822374 ,  0.5580728 ],
           [ 0.15772933,  0.48593197,  0.55801347],
           [ 0.15626973,  0.4896237 ,  0.557936  ],
           [ 0.15481488,  0.49331293,  0.55783967],
           [ 0.15336445,  0.49700003,  0.55772371],
           [ 0.1519182 ,  0.50068529,  0.55758733],
           [ 0.15047605,  0.50436904,  0.55742968],
           [ 0.14903918,  0.50805136,  0.5572505 ],
           [ 0.14760731,  0.51173263,  0.55704861],
           [ 0.14618026,  0.51541316,  0.55682271],
           [ 0.14475863,  0.51909319,  0.55657181],
           [ 0.14334327,  0.52277292,  0.55629491],
           [ 0.14193527,  0.52645254,  0.55599097],
           [ 0.14053599,  0.53013219,  0.55565893],
           [ 0.13914708,  0.53381201,  0.55529773],
           [ 0.13777048,  0.53749213,  0.55490625],
           [ 0.1364085 ,  0.54117264,  0.55448339],
           [ 0.13506561,  0.54485335,  0.55402906],
           [ 0.13374299,  0.54853458,  0.55354108],
           [ 0.13244401,  0.55221637,  0.55301828],
           [ 0.13117249,  0.55589872,  0.55245948],
           [ 0.1299327 ,  0.55958162,  0.55186354],
           [ 0.12872938,  0.56326503,  0.55122927],
           [ 0.12756771,  0.56694891,  0.55055551],
           [ 0.12645338,  0.57063316,  0.5498411 ],
           [ 0.12539383,  0.57431754,  0.54908564],
           [ 0.12439474,  0.57800205,  0.5482874 ],
           [ 0.12346281,  0.58168661,  0.54744498],
           [ 0.12260562,  0.58537105,  0.54655722],
           [ 0.12183122,  0.58905521,  0.54562298],
           [ 0.12114807,  0.59273889,  0.54464114],
           [ 0.12056501,  0.59642187,  0.54361058],
           [ 0.12009154,  0.60010387,  0.54253043],
           [ 0.11973756,  0.60378459,  0.54139999],
           [ 0.11951163,  0.60746388,  0.54021751],
           [ 0.11942341,  0.61114146,  0.53898192],
           [ 0.11948255,  0.61481702,  0.53769219],
           [ 0.11969858,  0.61849025,  0.53634733],
           [ 0.12008079,  0.62216081,  0.53494633],
           [ 0.12063824,  0.62582833,  0.53348834],
           [ 0.12137972,  0.62949242,  0.53197275],
           [ 0.12231244,  0.63315277,  0.53039808],
           [ 0.12344358,  0.63680899,  0.52876343],
           [ 0.12477953,  0.64046069,  0.52706792],
           [ 0.12632581,  0.64410744,  0.52531069],
           [ 0.12808703,  0.64774881,  0.52349092],
           [ 0.13006688,  0.65138436,  0.52160791],
           [ 0.13226797,  0.65501363,  0.51966086],
           [ 0.13469183,  0.65863619,  0.5176488 ],
           [ 0.13733921,  0.66225157,  0.51557101],
           [ 0.14020991,  0.66585927,  0.5134268 ],
           [ 0.14330291,  0.66945881,  0.51121549],
           [ 0.1466164 ,  0.67304968,  0.50893644],
           [ 0.15014782,  0.67663139,  0.5065889 ],
           [ 0.15389405,  0.68020343,  0.50417217],
           [ 0.15785146,  0.68376525,  0.50168574],
           [ 0.16201598,  0.68731632,  0.49912906],
           [ 0.1663832 ,  0.69085611,  0.49650163],
           [ 0.1709484 ,  0.69438405,  0.49380294],
           [ 0.17570671,  0.6978996 ,  0.49103252],
           [ 0.18065314,  0.70140222,  0.48818938],
           [ 0.18578266,  0.70489133,  0.48527326],
           [ 0.19109018,  0.70836635,  0.48228395],
           [ 0.19657063,  0.71182668,  0.47922108],
           [ 0.20221902,  0.71527175,  0.47608431],
           [ 0.20803045,  0.71870095,  0.4728733 ],
           [ 0.21400015,  0.72211371,  0.46958774],
           [ 0.22012381,  0.72550945,  0.46622638],
           [ 0.2263969 ,  0.72888753,  0.46278934],
           [ 0.23281498,  0.73224735,  0.45927675],
           [ 0.2393739 ,  0.73558828,  0.45568838],
           [ 0.24606968,  0.73890972,  0.45202405],
           [ 0.25289851,  0.74221104,  0.44828355],
           [ 0.25985676,  0.74549162,  0.44446673],
           [ 0.26694127,  0.74875084,  0.44057284],
           [ 0.27414922,  0.75198807,  0.4366009 ],
           [ 0.28147681,  0.75520266,  0.43255207],
           [ 0.28892102,  0.75839399,  0.42842626],
           [ 0.29647899,  0.76156142,  0.42422341],
           [ 0.30414796,  0.76470433,  0.41994346],
           [ 0.31192534,  0.76782207,  0.41558638],
           [ 0.3198086 ,  0.77091403,  0.41115215],
           [ 0.3277958 ,  0.77397953,  0.40664011],
           [ 0.33588539,  0.7770179 ,  0.40204917],
           [ 0.34407411,  0.78002855,  0.39738103],
           [ 0.35235985,  0.78301086,  0.39263579],
           [ 0.36074053,  0.78596419,  0.38781353],
           [ 0.3692142 ,  0.78888793,  0.38291438],
           [ 0.37777892,  0.79178146,  0.3779385 ],
           [ 0.38643282,  0.79464415,  0.37288606],
           [ 0.39517408,  0.79747541,  0.36775726],
           [ 0.40400101,  0.80027461,  0.36255223],
           [ 0.4129135 ,  0.80304099,  0.35726893],
           [ 0.42190813,  0.80577412,  0.35191009],
           [ 0.43098317,  0.80847343,  0.34647607],
           [ 0.44013691,  0.81113836,  0.3409673 ],
           [ 0.44936763,  0.81376835,  0.33538426],
           [ 0.45867362,  0.81636288,  0.32972749],
           [ 0.46805314,  0.81892143,  0.32399761],
           [ 0.47750446,  0.82144351,  0.31819529],
           [ 0.4870258 ,  0.82392862,  0.31232133],
           [ 0.49661536,  0.82637633,  0.30637661],
           [ 0.5062713 ,  0.82878621,  0.30036211],
           [ 0.51599182,  0.83115784,  0.29427888],
           [ 0.52577622,  0.83349064,  0.2881265 ],
           [ 0.5356211 ,  0.83578452,  0.28190832],
           [ 0.5455244 ,  0.83803918,  0.27562602],
           [ 0.55548397,  0.84025437,  0.26928147],
           [ 0.5654976 ,  0.8424299 ,  0.26287683],
           [ 0.57556297,  0.84456561,  0.25641457],
           [ 0.58567772,  0.84666139,  0.24989748],
           [ 0.59583934,  0.84871722,  0.24332878],
           [ 0.60604528,  0.8507331 ,  0.23671214],
           [ 0.61629283,  0.85270912,  0.23005179],
           [ 0.62657923,  0.85464543,  0.22335258],
           [ 0.63690157,  0.85654226,  0.21662012],
           [ 0.64725685,  0.85839991,  0.20986086],
           [ 0.65764197,  0.86021878,  0.20308229],
           [ 0.66805369,  0.86199932,  0.19629307],
           [ 0.67848868,  0.86374211,  0.18950326],
           [ 0.68894351,  0.86544779,  0.18272455],
           [ 0.69941463,  0.86711711,  0.17597055],
           [ 0.70989842,  0.86875092,  0.16925712],
           [ 0.72039115,  0.87035015,  0.16260273],
           [ 0.73088902,  0.87191584,  0.15602894],
           [ 0.74138803,  0.87344918,  0.14956101],
           [ 0.75188414,  0.87495143,  0.14322828],
           [ 0.76237342,  0.87642392,  0.13706449],
           [ 0.77285183,  0.87786808,  0.13110864],
           [ 0.78331535,  0.87928545,  0.12540538],
           [ 0.79375994,  0.88067763,  0.12000532],
           [ 0.80418159,  0.88204632,  0.11496505],
           [ 0.81457634,  0.88339329,  0.11034678],
           [ 0.82494028,  0.88472036,  0.10621724],
           [ 0.83526959,  0.88602943,  0.1026459 ],
           [ 0.84556056,  0.88732243,  0.09970219],
           [ 0.8558096 ,  0.88860134,  0.09745186],
           [ 0.86601325,  0.88986815,  0.09595277],
           [ 0.87616824,  0.89112487,  0.09525046],
           [ 0.88627146,  0.89237353,  0.09537439],
           [ 0.89632002,  0.89361614,  0.09633538],
           [ 0.90631121,  0.89485467,  0.09812496],
           [ 0.91624212,  0.89609127,  0.1007168 ],
           [ 0.92610579,  0.89732977,  0.10407067],
           [ 0.93590444,  0.8985704 ,  0.10813094],
           [ 0.94563626,  0.899815  ,  0.11283773],
           [ 0.95529972,  0.90106534,  0.11812832],
           [ 0.96489353,  0.90232311,  0.12394051],
           [ 0.97441665,  0.90358991,  0.13021494],
           [ 0.98386829,  0.90486726,  0.13689671],
           [ 0.99324789,  0.90615657,  0.1439362 ]]

    return LinearSegmentedColormap.from_list('viridis', cm_data)

viridis = load_viridis()
class Mesh():
    """
    Mesh class, for keeping all the mesh data.
    Input: directory [string] -> place where the domain files are.
    """
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            domain_x = np.loadtxt(directory+"domain_x.dat")
        except IOError:
            print "IOError with domain_x.dat"
        try:
            #We have to avoid ghost cells!
            domain_y = np.loadtxt(directory+"domain_y.dat")[2:-2]
        except IOError:
            print "IOError with domain_y.dat"
        self.xm = domain_x #X-Edge
        self.ym = domain_y #Y-Edge

        self.xmed = 0.5*(self.xm[:-1] + self.xm[1:]) #X-Center
        self.ymed = 0.5*(self.ym[:-1] + self.ym[1:]) #Y-Center

        self.dx = diff(self.xm)
        self.dy = diff(self.ym)
        self.surfy = outer(self.ym[:-1],self.dx)
        self.surfx = outer(self.dy,ones(self.dx.shape))

        self.vol = outer(.5*(self.ym[1:]**2-self.ym[:-1]**2),self.dx)

        self.nx = len(self.xmed)
        self.ny = len(self.ymed)-2

        #(Surfaces taken from the edges)
        #First we make 2D arrays for x & y, that are (theta,r)
#        T,R = meshgrid(self.xm, self.ym)
#        R2  = R*R
#        self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])


class Parameters():
    """
    Class for reading the simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, directory=''):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            params = open(directory+"variables.par",'r') #Opening the parfile
        except IOError:                  # Error checker.
            print  paramfile + " not found."
            return
        lines = params.readlines()     # Reading the parfile
        params.close()                 # Closing the parfile
        par = {}                       # Allocating a dictionary
        for line in lines:             #Iterating over the parfile
            name, value = line.split() #Spliting the name and the value (first blank)
            try:
                float(value)           # First trying with float
            except ValueError:         # If it is not float
                try:
                    int(value)         #                   we try with integer
                except ValueError:     # If it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # Filling the dictory
        self._params = par             # A control atribute, actually not used, good for debbuging
        for name in par:               # Iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) #Making the atributes at runtime


class Field(Mesh, Parameters):
    """
    Field class, it stores all the mesh, parameters and scalar data
    for a scalar field.
    Input: field [string] -> filename of the field
           staggered='c' [string] -> staggered direction of the field.
                                      Possible values: 'x', 'y', 'xy', 'yx'
           directory='' [string] -> where filename is
           dtype='float64' (numpy dtype) -> 'float64', 'float32',
                                             depends if FARGO_OPT+=-DFLOAT is activated
    """
    Q_dict = {'dens':'gasdens{0:d}.dat','vx':'gasvx{0:d}.dat','vy':'gasvy{0:d}.dat','momx':['dens','vx'],'momy':['dens','vy'],'pres':['dens']}
    name_dict={'dens':'$\\Sigma$', 'vx': '$v_\\phi$', 'vy': '$v_r$'}
    def __init__(self, Q, num, staggered='c', directory='', dtype='float64'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        Mesh.__init__(self, directory) #All the Mesh attributes inside Field!
        Parameters.__init__(self, directory) #All the Parameters attributes inside Field!
        self.staggered = staggered
        #Now, the staggering:
        if staggered.count('x')>0:
            self.x = self.xm[:-1] #Do not dump the last element
        else:
            self.x = self.xmed
        if staggered.count('y')>0:
            self.y = self.ym[:-1]
        else:
            self.y = self.ymed

        try:
            field = self.Q_dict[Q].format(num)
        except KeyError:
            print '%s not a valid Field'%Q
            print 'Please choose one of:, ',self.Q_dict.keys()
            raise

        self.name = Q
        self.math_name = self.name_dict[Q]
        self.data = self.__open_field(directory+field.format(num),dtype) #The scalar data is here.
        self.data = self.set_boundary()
        self.avg = self.data.mean(axis=1)
        self.wkz = self.ymin + (self.ymax-self.ymin)*self.wkzin
        self.wkzr = self.ymax - (self.ymax-self.ymin)*self.wkzout
        self.ft = fft.rfft(self.data,axis=1)/self.nx


    def __open_field(self, f, dtype):
        """
        Reading the data
        """
        try:
            field = fromfile(f, dtype=dtype)
        except IOError:
            print "Couldn't find %s, trying again with un-merged file"%f
            try:
                field = fromfile(f.replace('.dat','_0.dat'),dtype=dtype)
            except IOError:
                raise

        return field.reshape(self.ny, self.nx)


    def recalculate(self):
        self.avg = self.data.mean(axis=1)
        self.ft = fft.rfft(self.data,axis=1)/self.nx

    def grad(self):
        q = copy.copy(self.data)
        res = zeros(q.shape)
        one_dim = False
        try:
            q.shape[1]
        except IndexError:
            one_dim = True


        for i in range(1,q.shape[0]-1):
            if one_dim:
                res[i] = (q[i+1]-q[i-1])/(self.y[i+1]-self.y[i-1])
            else:
                res[i,:] = (q[i+1,:]-q[i-1,:])/(self.y[i+1]-self.y[i-1])

        if one_dim:
            res[0] = (q[1]-q[0])/(self.y[1]-self.y[0])
            res[-1] = (q[-1]-q[-2])/(self.y[-1]-self.y[-2])
        else:
            res[0,:] = (q[1,:]-q[0,:])/(self.y[1]-self.y[0])
            res[-1,:] = (q[-1,:]-q[-2,:])/(self.y[-1]-self.y[-2])
        return res


    def set_boundary(self,accretion=True):
        if self.name == 'vx':
            iact = 1; igh = 0;
            inner_ring = (self.data[iact,:] + self.ymed[iact]*self.omegaframe)*sqrt(self.ymed[iact]/self.ymed[igh]) - self.ymed[igh]*self.omegaframe
            iact = -2; igh = -1;
            outer_ring = (self.data[iact,:] + self.ymed[iact]*self.omegaframe)*sqrt(self.ymed[iact]/self.ymed[igh]) - self.ymed[igh]*self.omegaframe
        if self.name == 'vy':
            if accretion:
                iact=1; igh=0;
                inner_ring = -1.5*self.alpha*self.aspectratio*self.aspectratio*self.ymed[igh]**(2*self.flaringindex-.5)
                inner_ring *= ones(self.data[igh,:].shape)
                iact=-2; igh=-1;
                outer_ring = self.data[iact,:]
            else:
                inner_ring = self.data[iact,:]
                outer_ring = self.data[iact,:]

        if self.name == 'dens':
            if accretion:
                iact=1;igh=0;
                inner_ring = self.data[iact,:]* (self.ymed[iact]/self.ymed[igh])**(2*self.flaringindex+.5)
                iact=-2;igh=-1;
                nu_0 = self.alpha*self.aspectratio*self.aspectratio
                fac = 3*pi*nu_0* (self.ymed[iact])**(2*self.flaringindex+1)
                outer_ring = (self.data[iact,:]*fac  + self.mdot*(sqrt(self.ymed[igh])-sqrt(self.ymed[iact])))/(3*pi*nu_0*self.ymed[igh]**(2*self.flaringindex+1))
            else:
                iact=1;igh=0;
                inner_ring = self.data[iact,:]
                iact=-2;igh=-1;
                outer_ring = self.data[iact,:]

        return vstack( (inner_ring, self.data, outer_ring) )

    def center_to_edge(self,direction='y'):
        newfield = copy.deepcopy(self)
        if direction.count('y') > 0:
            newfield.data = vstack( (newfield.data, newfield.data[-1,:]) )
            newfield.data = .5*(newfield.data[1:,:]+  newfield.data[:-1,:])
            newfield.staggered='y'
        if direction.count('x') > 0:
            newfield.data = hstack( (newfield.data[:,-1].reshape(len(newfield.data[:,-1]),1), newfield.data, newfield.data[:,0].reshape(len(newfield.data[:,0]),1)) )
            newfield.data = .5*(newfield.data[:,1:]+  newfield.data[:,:-1])
            newfield.staggered += 'x'
        newfield.avg = newfield.data.mean(axis=1)
        return newfield
    def edge_to_center(self):
        newfield = copy.deepcopy(self)
        if newfield.staggered.count('y') > 0:
            newfield.data = vstack( (newfield.data, newfield.data[-1,:]) )
            newfield.data = .5*(newfield.data[1:,:] + newfield.data[:-1,:])
            newfield.staggered.strip('y')
        if newfield.staggered.count('x') > 0:
            newfield.data = hstack( (newfield.data[:,-1].reshape(len(newfield.data[:,-1]),1), newfield.data, newfield.data[:,0].reshape(len(newfield.data[:,0]),1)) )
            newfield.data = .5*(newfield.data[:,1:] + newfield.data[:,:-1])
            newfield.staggered.strip('x')
        newfield.avg = newfield.data.mean(axis=1)
        return newfield


#    def center_to_edge(self,direction='y'):
#        dqm = zeros(rho.shape)
#        dqp = zeros(rho.shape)
#        slope = zeros(rho.shape)
#        mdot = zeros(rho.shape)
#
#
#        dqm[1:-1,:] =(rho[1:-1,:]  - rho[:-2,:])/dy[1:-1,:]
#        dqp[1:-1,:] =(rho[2:,:]  - rho[1:-1,:])/dy[2:,:]
#        ind = sign(dqm*dqp)
#        slope = (ind>0).astype(int) * 2*dqm*dqp/(dqm+dqp)
#        mdot[1:-1,:] = (vy>0).astype(int)*(rho[:-2,:]+.5*slope[:-2,:]*dy[:-2,:]) + (vy<=0).astype(int)*(rho[1:,:]-.5*slope[1:,:]*dy[1:,:])

    def plotmode(self,m,ax=None,norm=1,shift=0,planet=None,xlims=None,ylims=None,logx=False,**karg):

        try:
            m[0]
        except TypeError:
            m= array([m])


        if planet is None:
            y = copy.copy(self.y)
            xstr = '$r$'
        else:
            y = (self.y-planet)/(self.aspectratio*planet)
            xstr = '$(r-a)/H$'

        fontsize=karg.pop('fontsize',20)
        figsize = karg.pop('figsize',(10,8))
        if ax is None:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        color = karg.pop('color','k')
        for i in m:
            ax.plot(y,self.ft[:,i].real,'-',label='$m=%d$'%i,**karg)
        ax.set_prop_cycle(None)
        for i in m:
            ax.plot(y,self.ft[:,i].imag,'--',**karg)

        ax.set_xlabel(xstr,fontsize=fontsize)
        ax.set_ylabel('$\\hat{' + self.math_name.strip('$') + '}$' ,fontsize=fontsize)
        if len(m) < 6:
            ax.set_title('$m=$' + ','.join(['%d'%i for i in m]),fontsize=fontsize)
            ax.legend(loc='best')
        else:
            ax.set_title('$m = %d...%d$'%(min(m),max(m)),fontsize=fontsize)

        if logx and planet is None:
            ax.set_xscale('log')

        ax.minorticks_on()
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

    def plotavg(self,ax=None,log=False,logx=False,logy=False,xlims=None,ylims=None,planet=None,norm=1,shift=0,**karg):
        fontsize=karg.pop('fontsize',20)
        figsize = karg.pop('figsize',(10,8))


        if ax is None:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        if self.staggered.count('y')>0:
            y = copy.copy(self.ym[:-1])
        else:
            y = copy.copy(self.ymed)

        xstr = '$r$'
        if planet is not None:
            y = (y - planet)/(self.aspectratio*planet)
            xstr = '$ (r-a)/H $'
            logx=False
            if log:
                logy=True
                log=False

        ax.plot(y,(self.avg-shift)/norm,**karg)
        if logx or log:
            ax.set_xscale('log'),
        if logy or log:
            ax.set_yscale('log')
        ax.set_xlabel(xstr,fontsize=fontsize)
        ax.set_ylabel(self.math_name,fontsize=fontsize)
        ax.minorticks_on()
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)



    def plot(self, norm=1,shift=0,ax=None,log=False,logx=False, abslog=False,cartesian=False, cmap=viridis, **karg):
        """
        A layer to plt.imshow or pcolormesh function.
        if cartesian = True, pcolormesh is launched.
        """

#        if self.x[0] == pi or self.x[-1] != pi:
#            dp = diff(self.x)[0]
#            dp *= 2
#            yn = self.data[:,-1] + (pi - self.x[-1])*(self.data[:,0]-self.data[:,-1])/dp
#            data = copy.copy(self.data)
#            data[:,-1] = yn
#            data[:,0] = yn
#            x = copy.copy(self.x)
#            x[-1] = pi
#            x[0] = -pi
#
#        x = copy(self.x)
        data = (copy.copy(self.data)-shift)/norm
        x =copy.copy(self.x)
        shading = karg.pop('shading','gouraud')
        interpolation = karg.pop('interpolation','bilinear')
        fontsize = karg.pop('fontsize',20)
        xlims = karg.pop('xlims',None)
        ylims = karg.pop('ylims',None)
        T,R = meshgrid(x,self.y)
        if ax == None:
            fig=figure()
            ax=fig.add_subplot(111)
        if log:
            data = np.log10(data)
        if abslog:
            data = np.log10(np.abs(data))
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            line2d=ax.pcolormesh(X,Y,data,cmap=cmap,shading=shading,**karg)
            if xlims is not None:
                ax.set_xlim(xlims)
            if ylims is not None:
                ax.set_ylim(ylims)

        else:
            if logx:
                line2d=ax.pcolormesh(log10(R),T,data, cmap = cmap,shading=shading,**karg)
                ax.set_xlabel('$\\log_{10}(r)$',fontsize=fontsize)
                if xlims is not None:
                    ax.set_xlim(log10(xlims[0]),log10(xlims[1]))
                else:
                    ax.set_xlim(log10(self.y[0]),log10(self.y[-1]))
            else:
                line2d=ax.pcolormesh(R,T,data, cmap = cmap,shading=shading,**karg)
                ax.set_xlabel('$r$',fontsize=fontsize)
                if xlims is not None:
                    ax.set_xlim(xlims)
                else:
                    ax.set_xlim(self.y[0],self.y[-1])

            #origin='lower',aspect='auto',interpolation=interpolation,extent=[x[0],x[-1],self.y[0],self.y[-1]],**karg)
            ax.set_ylabel('$\\phi$',fontsize=fontsize)
            if ylims is not None:
                ax.set_ylim(ylims)
            else:
                ax.set_ylim(self.x[0],self.x[-1])
        cbar = colorbar(line2d,ax=ax)
        if log:
            ax.set_title('$\\log_{10}$'+self.math_name,fontsize=fontsize)
        elif abslog:
            ax.set_title('$\\log_{10}|$'+self.math_name+'$|$',fontsize=fontsize)
        else:
            ax.set_title(self.math_name,fontsize=fontsize)

    def contour(self, log=False, cartesian=False, **karg):
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        ax = gca()
        T,R = meshgrid(self.x,self.y)
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            contour(X,Y,data,**karg)
        else:
            contour(T,R,data,**karg)


    def shift_field(Field,direction):
        import copy
        """
        Half cell shifting along the direction provided by direction
        direction can be ('x','y', 'xy', 'yx').

        After a call of this function, Field.xm/xmed has not
        sense anymore (it is not hard to improve).
        """
        F = copy.deepcopy(Field)
        if direction.count('x')>0:
            F.data = 0.5*(Field.data[:,1:]+Field.data[:,:-1])
            F.x = 0.5*(Field.x[1:]+Field.x[:-1])
        if direction.count('y')>0:
            F.data = 0.5*(F.data[1:,:]+F.data[:-1,:])
            F.y = 0.5*(F.y[1:]+F.y[:-1])

        F.nx = len(F.x)
        F.ny = len(F.y)

        return F


    def cut_field(Field, direction, side):
        """
        Cutting a field:
        Input: field --> a Field class
               axis  --> 'x', 'y' or 'xy'
               side  --> 'p' (plus), 'm' (minnus), 'pm' (plus/minnus)
        """

        cutted_field = copy.deepcopy(Field)
        ny,nx = Field.ny, Field.nx
        mx = my = px = py = 0

        if direction.count('x')>0:
            if side.count('m')>0:
                mx = 1
            if side.count('p')>0:
                px = 1
        if direction.count('y')>0:
            if side.count('m')>0:
                my = 1
            if side.count('p')>0:
                py = 1

        cutted_field.data = Field.data[my:ny-py,mx:nx-px]
        cutted_field.x = cutted_field.x[mx:nx-px]
        cutted_field.y = cutted_field.y[my:ny-py]

        return cutted_field

    def vector_field(vx,vy, **karg):
        nsx = nsy = 3
        T,R = meshgrid(vx.x[::nsx],vx.y[::nsy])
        X = R*cos(T)
        Y = R*sin(T)
        vx = vx.data[::nsy,::nsx]
        vy = vy.data[::nsy,::nsx]
        U = vy*cos(T) - vx*sin(T)
        V = vy*sin(T) + vx*cos(T)
        ax = gca()
        ax.quiver(X,Y,U,V,scale=5,pivot='midle', **karg)

    def euler(vx, vy, x, y, reverse):
        """
        Euler integrator for computing the streamlines.
        Parameters:
        ----------

        x,y: Float.
             starter coordinates.
        reverse: Boolean.
                 If reverse is true, the integratin step is negative.
        Reverse inverts the sign of the velocity

        Output
        ------

        (dx,dy): (float,float).
                 Are the azimutal and radial increment.
                 Only works for cylindrical coordinates.
        """
        sign = 1.0
        if reverse:
            sign = -1
        vphi = get_v(vx,x,y)
        vrad = get_v(vy,x,y)
        if vphi == None or vrad == None: #Avoiding problems...
            return None,None
        l = np.min((((vx.xmax-vx.xmin)/vx.nx),((vx.ymax-vx.ymin)/vx.ny)))
        h = 0.5*l/np.sqrt((vphi**2+vrad**2))

#        return sign*h*np.array([vphi/y,vrad])
        return sign*h*np.array([vphi,vrad])


    def bilinear(x,y,f,p):
        """
        Computing bilinear interpolation.
        Parameters
        ----------
        x = (x1,x2); y = (y1,y2)
        f = (f11,f12,f21,f22)
        p = (x,y)
        where x,y is the interpolated point and
        fij is the value of the function at the
        point (xi,yj).
        Output
        ------
        f(p): Float.
              The interpolated value of the function f(p) = f(x,y)
        """

        xp  = p[0]; yp   = p[1]; x1  = x[0]; x2  = x[1]
        y1  = y[0]; y2  = y[1];  f11 = f[0]; f12 = f[1]
        f21 = f[2]; f22 = f[3]
        t = (xp-x1)/(x2-x1);    u = (yp-y1)/(y2-y1)

        return (1.0-t)*(1.0-u)*f11 + t*(1.0-u)*f12 + t*u*f22 + u*(1-t)*f21

    def get_v(v, x, y):
        """
        For a real set of coordinates (x,y), returns the bilinear
        interpolated value of a Field class.
        """

        i = int((x-v.xmin)/(v.xmax-v.xmin)*v.nx)
        j = int((y-log(v.ymin))/(log(v.ymax/v.ymin))*v.ny)

        if i<0 or j<0 or i>v.data.shape[1]-2 or j>v.data.shape[0]-2:
            return None

        f11 = v.data[j,i]
        f12 = v.data[j,i+1]
        f21 = v.data[j+1,i]
        f22 = v.data[j+1,i+1]
        try:
            x1  = v.x[i]
            x2  = v.x[i+1]
            y1  = log(v.y[j])
            y2  = log(v.y[j+1])
            return bilinear((x1,x2),(y1,y2),(f11,f12,f21,f22),(x,y))
        except IndexError:
            return None



    def get_stream(vx, vy, x0, y0, nmax=1000000, maxlength=4*np.pi, bidirectional=True, reverse=False):
        """
        Function for computing a streamline.
        Parameters:
        -----------

        x0,y0: Floats.
              Initial position for the stream
        nmax: Integer.
              Maxium number of iterations for the stream.
        maxlength: Float
                   Maxium allowed length for a stream
        bidirectional=True
                      If it's True, the stream will be forward and backward computed.
        reverse=False
                The sign of the stream. You can change it mannualy for a single stream,
                but in practice, it's recommeneded to use this function without set reverse
                and setting bidirectional = True.

        Output:
        -------

        If bidirectional is False, the function returns a single array, containing the streamline:
        The format is:

                                          np.array([[x],[y]])

        If bidirectional is True, the function returns a tuple of two arrays, each one with the same
        format as bidirectional=False.
        The format in this case is:

                                (np.array([[x],[y]]),np.array([[x],[y]]))

        This format is a little bit more complicated, and the best way to manipulate it is with iterators.
        For example, if you want to plot the streams computed with bidirectional=True, you should write
        some similar to:

        stream = get_stream(x0,y0)
        ax.plot(stream[0][0],stream[0][1]) #Forward
        ax.plot(stream[1][0],stream[1][1]) #Backward

        """

        if bidirectional:
            s0 = get_stream(vx, vy, x0, y0, reverse=False, bidirectional=False, nmax=nmax,maxlength=maxlength)
            s1 = get_stream(vx, vy, x0, y0, reverse=True,  bidirectional=False, nmax=nmax,maxlength=maxlength)
            return (s0,s1)

        l = 0
        x = [x0]
        y = [y0]

        for i in xrange(nmax):
            ds = euler(vx, vy, x0, y0, reverse=reverse)
            if ds[0] == None:
                if(len(x)==1):
                    print "There was an error getting the stream, ds was NULL (see get_stream)."
                break
            l += np.sqrt(ds[0]**2+ds[1]**2)
            dx = ds[0]
            dy = ds[1]
            if(np.sqrt(dx**2+dy**2)<1e-13):
                print "Warning (get_stream): ds is very small, maybe you're in a stagnation point."
                print "Try selecting another initial point."
                break
            if l > maxlength:
                print "maxlength reached: ", l
                break
            x0 += dx
            y0 += dy
            x.append(x0)
            y.append(y0)
        return np.array([x,y])

    def get_random_streams(vx, vy, xmin=None, xmax=None, ymin=None, ymax=None, n=30, nmax=100000):
        if xmin == None:
            xmin = vx.xmin
        if ymin == None:
            ymin = vx.ymin
        if xmax == None:
            xmax = vx.xmax
        if ymax == None:
            ymax = vx.ymax
        X = xmin + np.random.rand(n)*(xmax-xmin)
        Y = ymin + np.random.rand(n)*(ymax-ymin)
        streams = []
        counter = 0
        for x,y in zip(X,Y):
            stream = get_stream(vx, vy, x, y, nmax=nmax, bidirectional=True)
            streams.append(stream)
            counter += 1
            print "stream ",counter, "done"
        return streams

    def plot_random_streams(streams, cartesian=False, **kargs):
        ax = plt.gca()
        print np.shape(streams)
        for stream in streams:
            for sub_stream in stream:
                if cartesian:
                    ax.plot(sub_stream[1]*cos(sub_stream[0]),sub_stream[1]*sin(sub_stream[0]),**kargs)
                else:
                    ax.plot(sub_stream[0],sub_stream[1],**kargs)

from scipy.interpolate import interp1d,interp2d


class Sim(Mesh,Parameters):
    def __init__(self,i,directory='',p=0):
        if directory != '':
            if directory[-1] != '/':
                directory += '/'
        self.directory = directory
        self.dens = Field('dens',i,directory=directory)
        self.vp = Field('vx',i,directory=directory,staggered='x')
        self.vr = Field('vy',i,directory=directory,staggered='y')

        Mesh.__init__(self,directory)


        self.mdot,self.rhos,self.rhosp = self.calc_flux()
        self.vrf = self.mdot.avg/(-2*pi*self.mdot.y*self.rhos.avg)
        self.sig0  = self.dens.mdot/(3*pi*self.nu(self.dens.y))
#        self.vp = self.vp.shift_field('x')
#        self.vr = self.vr.shift_field('y')
#        self.vp = self.vp.cut_field(direction='y',side='p')
#        self.vr = self.vr.cut_field(direction='x',side='p')
#        self.dens = self.dens.cut_field(direction='xy', side='p')
        self.r  = self.dens.y
        self.phi = self.dens.x
        self.dphi = self.dx
#        self.dlr = diff(log(self.r))[0]
#        self.dr = self.dlr*self.r
        self.dr = self.dy
        self.dlr = self.dy/self.ymed

        self.dlr = self.dlr[:self.dens.ny-1]
        self.dr = self.dr[:self.dens.ny-1]



#        self.pp, self.rr = meshgrid(self.phi,self.r)
#    	self.dbar = self.dens.data.mean(axis=1)
#        self.vrbar0 = self.vr.data.mean(axis=1)
#    	self.vrbar = (self.dens.data*self.vr.data).mean(axis=1)
#    	self.vpbar = (self.dens.data*self.vp.data).mean(axis=1)
#    	self.mdot_full = -2*pi*self.rr*(self.vr.data*self.dens.data)
#        self.mdot = self.mdot_full.mean(axis=1)
        try:
    	    _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,self.t,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))[i,:]
        except IndexError:
            try:
    	        _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,self.t,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))[-1,:]
            except IndexError:
    	        _,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,self.t,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))


    	self.a = sqrt(self.px**2  + self.py**2)
        self.K = self.mp**2 / (self.dens.alpha * self.dens.aspectratio**5)

        self.vp.data += self.vp.y[:,newaxis]*self.omf
        self.vp.recalculate()

        try:
            self.nu0 = self.dens.alpha*self.dens.aspectratio*self.dens.aspectratio
        except AttributeError:
            self.nu0 = self.dens.alphaviscosity*self.dens.aspectratio*self.dens.aspectratio


    	self.tvisc = self.dens.ymax**2/self.nu(self.dens.ymax)
        self.tviscp = self.a**2/(self.nu(self.a))
        self.torb = 2*np.pi * self.a**(1.5)
        self.rh = (self.mp/3)**(1./3) * self.a

        self.dTr,self.Lambda,self.Fh=self.calc_torques()

##    	self.omega = zeros(self.vp.data.shape)
##    	self.vpc = zeros(self.vp.data.shape)
##        self.l = zeros(self.vp.data.shape)
##        for i in range(len(self.r)):
##            self.omega[i,:] = self.vp.data[i,:]/self.r[i] + self.a**(-1.5)
##            self.vpc[i,:] = self.omega[i,:]*self.r[i]
##            self.l[i,:] = self.r[i] * self.vpc[i,:]
##
##        self.Fj_full = -2*pi*self.dens.data*self.nu(self.rr)*self.rr**3 * self.grad(self.omega)
##        self.Fj = self.Fj_full.mean(axis=1)
##        ind_excl = (self.r>=self.a*(1-self.dens.thicknesssmoothing*self.dens.aspectratio))&(self.r<=self.a*(1+self.dens.thicknesssmoothing*self.dens.aspectratio))
##        self.dTr = (-self.rr*self.dens.data * self.dp_potential(self.rr,self.pp))
##        self.dTr_excl = self.dTr.mean(axis=1)
##        self.dTr_excl[ind_excl] = 0
##        self.dTr_tot_excl = 2*pi*trapz(self.dTr_excl* self.dr)
##        self.dTr_mean = self.dTr.mean(axis=1)
##        self.dTr_total = 2*pi*trapz(self.dTr_mean*self.dr)
##        self.dTr_out = self.dphi*trapz((self.dTr.sum(axis=1)*self.dr)[self.r>=self.a],x=self.r[self.r>=self.a])
##        self.dTr_in = self.dphi*trapz((self.dTr.sum(axis=1)*self.dr)[self.r<=self.a],x=self.r[self.r<=self.a])
##
##        self.mdotdlr = (self.mdot_full * self.grad(self.l)).mean(axis=1)
##
##        self.safety_fac  = .5
##        try:
##            self.dbar0 = self.dens.mdot/(3*pi*self.nu(self.r))
##        except AttributeError:
##            pass
##        self.vr_visc = -1.5*self.nu(self.r)/self.r
##        self.steady_state_calc()
##        self.mdoti,self.dbar0=self.inner_disk_sol((self.dens.ymin*1.2,self.a*.5))
##
##        self.A = self.dTr_total/self.mdoti
##        self.A_excl = self.dTr_tot_excl/self.mdoti
    def fpred(self,x):
        scalar = False
        try:
            x[1]
        except IndexError:
            scalar = True
        ind = self.r<=x
        res = self.dTr_mean*2*pi + self.mdotdlr
        if scalar:
            return trapz( (res*self.dr)[ind])
        else:
            return array([trapz( (res*self.dr)[self.r<=i]) for i in x])
    def calc_torques(self):
        res_dTr = copy.deepcopy(self.rhos)

        xx,yy = meshgrid(res_dTr.x,res_dTr.y)
        res_dTr.data *=  (-yy * self.dp_potential(yy,xx))

        res_dTr.name = 'Lambda_ex'
        res_dTr.math_name = '$\\Lambda_{ex}$'
        res_dTr.recalculate()


        res_fH = copy.deepcopy(self.rhosp)
        res_dep = copy.deepcopy(self.dens)
        # Assume a keplerian rotation for now
        res_fH.data *=  (3*pi*self.nu(yy) * yy**(.5))
        res_fH.name = 'Fh'
        res_fH.math_name = '$F_{H,\\nu}$'
        res_fH.recalculate()
        res_dep.data  = -self.mdot.data /(2*sqrt(yy)) + res_fH.grad()
        res_dep.name = 'Lambda_dep'
        res_dep.math_name = '$\\Lambda_{dep}$'
        res_dep.recalculate()

        return res_dTr,res_dep,res_fH

    def calc_flux(self):
        res = copy.deepcopy(self.vr)
        res1 = copy.deepcopy(self.vr)
        rho = copy.copy(self.dens.data)
        vy = copy.copy(self.vr.data)
        surfy = copy.copy(self.surfy)

        dqm = zeros(rho.shape)
        dqp = zeros(rho.shape)
        slope = zeros(rho.shape)
        mdot = zeros(rho.shape)
        rhos = copy.copy(self.dens.data)
        _,dyy = meshgrid(self.dx,self.dy)

        dqm[1:-1,:] = (rho[1:-1,:]-rho[:-2,:])/dyy[1:-1,:]
        dqp[1:-1,:] = (rho[2:,:]-rho[1:-1,:])/dyy[2:,:]
        slope[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
        rhos[1:-1,:] = (vy[1:-1,:]>0).astype(int) * (rho[:-2,:]+.5*slope[:-2,:]*dyy[:-2,:])
        rhos[1:-1,:] += (vy[1:-1,:]<=0).astype(int) * (rho[1:-1,:]-.5*slope[1:-1,:]*dyy[1:-1,:])
        mdot[1:-1,:] = -self.nx*vy[1:-1,:] * rhos[1:-1,:] * surfy[1:-1,:]
        mdot[0,:] = mdot[1,:]
        mdot[-1,:] = mdot[-2,:]

        res.data = mdot
        res.name = 'Mdot'
        res.math_name = '$\\dot{M}$'
        res1.data = rhos
        res1.name = 'Rhostar'
        res1.math_name = '$\\Sigma^*$'
        res.recalculate()
        res1.recalculate()

        # phi direction
        res2 = copy.deepcopy(self.vp)
        rho = copy.copy(self.dens.data)
        ns = (len(rho[:,0]),1)
        dx = self.dx[0]


        nrho = hstack( (rho[:,-1].reshape(ns),rho,rho[:,0].reshape(ns)) )

        dqm = zeros(nrho.shape)
        dqp = zeros(nrho.shape)
        slopep = zeros(nrho.shape)

        dqm[:,1:] = (nrho[:,1:]-nrho[:,:-1])/dx
        dqp[:,:-1] = (nrho[:,1:]-nrho[:,:-1])/dx
        dqm[:,0] = (nrho[:,0]-rho[:,-2])/dx
        dqp[:,-1] = (rho[:,1]-nrho[:,-1])/dx

        slopep = ( (dqm*dqp) > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))

        res2.data = nrho[:,:-2] + .5*slopep[:,:-2]*dx
#        res2.data[:,:-1] = .5*(self.dens.data[:,1:] + self.dens.data[:,:-1])
#        res2.data[:,-1] = .5*(self.dens.data[:,-1] + self.dens.data[:,0])
        res2.name = 'Rhostarp'
        res2.math_name = '$\\Sigma_\\phi^*$'
        res2.recalculate()


        return res,res1,res2
    def grad(self,q):
        res = zeros(q.shape)
        one_dim = False
        try:
            q.shape[1]
        except IndexError:
            one_dim = True


        for i in range(1,q.shape[0]-1):
            if one_dim:
                res[i] = (q[i+1]-q[i-1])/(self.r[i+1]-self.r[i-1])
            else:
                res[i,:] = (q[i+1,:]-q[i-1,:])/(self.r[i+1]-self.r[i-1])

        if one_dim:
            res[0] = (q[1]-q[0])/(self.r[1]-self.r[0])
            res[-1] = (q[-1]-q[-2])/(self.r[-1]-self.r[-2])
        else:
            res[0,:] = (q[1,:]-q[0,:])/(self.r[1]-self.r[0])
            res[-1,:] = (q[-1,:]-q[-2,:])/(self.r[-1]-self.r[-2])
        return res
    def nu(self,x):
        return self.dens.aspectratio**2 * self.dens.alpha  * x**(2*self.dens.flaringindex+.5)
    def vr_nu(self,x):
        return -1.5*self.nu(x)/x
    def scaleH(self,x):
        return self.dens.aspectratio * x**(self.dens.flaringindex + 1)

    def inner_disk_sol(self,rlims):
        ind = (self.r>=rlims[0])&(self.r<=rlims[1])
        popt,pcov = curve_fit(lambda x,a: a/(3*pi*self.nu(x)),self.r[ind],self.dbar[ind])
        return popt[0],popt[0]/(3*pi*self.nu(self.r))

    def plot2d(self,q,axex=None,fig=None,logscale=False,logr=True,rlims=None,plims=None,norm=False,**kargs):
        if fig == None:
            #fig,axes=subplots(1,2)
            figsize = kargs.pop('figsize',(20,15))
            fig = figure(figsize=figsize);
            gs = gridspec.GridSpec(3,4)
            axcart = fig.add_subplot(gs[:2,:-2])
            axcyl = fig.add_subplot(gs[:2,-2:])
            axdbar = fig.add_subplot(gs[-1,:])


        if q=='dens':
            fld=self.dens
            dat = copy.copy(self.dens.data.transpose())
            dstr0 = '$\\Sigma(r)$'
            dbar = copy.copy(self.dbar)
            if norm:
                dat0 = copy.copy(self.dbar0)
        elif q=='vr':
            fld =self.vr
            dat = copy.copy(self.vr.data.transpose())
            dstr0 = '$v_r$'

            dbar = copy.copy(self.vrbar)
            if norm:
                dat0 = copy.copy(self.vr_visc)
        elif q=='vp':
            fld=self.vp
            dat = copy.copy(self.vp.data.transpose())
            dbar = copy.copy(self.vpbar)
            dstr0 = '$v_\\phi$'
            if norm:
                dat0 = pow(self.r,-1.5)

        else:
            print '%s is not a valid option' % q
            return

        if logscale:
            dstr = '$\\log_{10}$'+ dstr0
        else:
            dstr = dstr0
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)

        print rlims,plims

        rinds = (self.r<=rlims[1])&(self.r>=rlims[0])
        pinds = (self.phi<=plims[1])&(self.phi>=plims[0])

        print self.phi[pinds]
        dat = dat[:,rinds][pinds,:]

        if norm:
            if q=='vp':
                dbar = (dbar- dat0)[rinds]
                for i in range(dat.shape[0]):
                    dat[i,:] -= dat0[rinds]
            else:
                dbar =  (dbar/dat0)[rinds]
                for i in range(dat.shape[0]):
                    dat[i,:] /= dat0[rinds]


        if logscale:
            dat = log10(dat)
        r = self.r[rinds]
        lr = log(r)
        phi = self.phi[pinds]

        print r,phi
        rlims0 = copy.copy(rlims)
        if logr:
            rstr = '$\\ln r$'
            rlims = (log(rlims[0]),log(rlims[1]))
            rr,pp = meshgrid(lr,phi)
        else:
            rstr = '$r$'
            rr,pp = meshgrid(r,phi)

        print rlims

        cmap = kargs.pop('cmap',viridis)
        line2d = axcyl.pcolormesh(rr,pp,dat,cmap=cmap,shading='gouraud')
        cbar = colorbar(line2d,ax=axcyl)
        cbar.set_label(dstr,fontsize=20)
        axcyl.set_xlim(rlims)
        axcyl.set_ylim(plims)
        axcyl.set_xlabel(rstr,fontsize=20)
        axcyl.set_ylabel('$\\phi$',fontsize=20)

        axdbar.plot(r,dbar)
        axdbar.set_xlabel('$r$',fontsize=20)
        axdbar.set_ylabel(dstr0,fontsize=20)
        if norm:
            axdbar.axhline(1,color='k')
        if logr:
            axdbar.set_xscale('log')
        if logscale:
            axdbar.set_yscale('log')
        axdbar.set_xlim(rlims0)
        fld.plot(cartesian=True,ax=axcart,log=logscale)


        return fig,[axcart,axcyl,axdbar]

    def potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + smoothing + self.a**2 - 2*r*self.a*cos(phi)
        rad = sqrt(rad)
        return -self.mp/rad
    def dp_potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + self.a**2 - 2*r*self.a*cos(phi) + smoothing
        rad = rad**(1.5)
        return self.mp * r*self.a*sin(phi)/rad
    def dr_potential(self,r,phi):
        if not self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(self.a)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + self.a**2 - 2*r*self.a*cos(phi) + smoothing
        rad = rad**(1.5)
        return self.mp * (r-self.a*cos(phi))/rad


    def conv_fargo_tq_1d(self,r,tq):
        tq *= -self.mp/(r**2 * self.dlr*2*pi)
        return tq
    def conv_fargo_tq_tot(self,tq):
        tq *= -self.mp/(2*pi)
        return tq

    def animate_mdot(self,irange,logspacing=True,cmap=viridis,**kargs):

        self.__init__(0,directory=self.directory)
        temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
        dens_list = [log10(temp0)]
        divm_list = [temp1]
        divm_mean_list = [temp2]
        drmx_mean_list = [temp3]
        dpmy_mean_list = [temp4]
        mx_mean_list = [temp5]
        my_mean_list = [temp6]
        t = [self.t]

        fig,axes = subplots(1,4,figsize=(15,10))

        axes[3].set_title('t = %.1f' % t[0],color='w')
        line0=axes[0].imshow(dens_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs)
        line1=axes[1].imshow(divm_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs);

        line2,=axes[2].plot(self.r,divm_mean_list[0],label='< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >')
        line3,=axes[2].plot(self.r,drmx_mean_list[0],label='<$\\nabla_r (r \\Sigma v_r)$>')
        line4,=axes[2].plot(self.r,dpmy_mean_list[0],label='-<$\\nabla_\\phi(\\Sigma v_\\phi)$>')
        axes[2].legend(loc='best')
        axes[2].set_ylim(-.0001,.0001)
        line5,=axes[3].plot(self.r,mx_mean_list[0],label='<-$r \\Sigma v_r$>')
        line6,=axes[3].plot(self.r,my_mean_list[0],label='-<$\\Sigma v_\\phi$>')
        axes[3].legend(loc='best')

        for i in irange:
            if i > 0:
                self.__init__(i,directory=self.directory)
                t.append(self.t)
                temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
                dens_list.append(log10(temp0).transpose())
                divm_list.append(temp1.transpose())
                divm_mean_list.append(temp2)
                drmx_mean_list.append(temp3)
                dpmy_mean_list.append(temp4)
                mx_mean_list.append(temp5)
                my_mean_list.append(temp6)


        for i in range(len(t)):
            axes[2].set_title('t = %.1f' % t[i])
            line0.set_data(dens_list[i])
            line1.set_data(divm_list[i])
            line2.set_ydata(divm_mean_list[i])
            line3.set_ydata(drmx_mean_list[i])
            line4.set_ydata(dpmy_mean_list[i])
            line5.set_ydata(mx_mean_list[i])
            line6.set_ydata(my_mean_list[i])
            fig.canvas.draw()


    def plot_mdot(self,logspacing=True):
        dens,divm, divm_mean, drmx_mean, dpmy_mean, mx_mean,my_mean = self.calculate_total_mdot(logspacing)
        fig,axes = subplots(1,4,figsize=(15,10))

        axes[0].imshow(log10(dens),aspect='auto',origin='lower');
        axes[1].imshow(divm,aspect='auto',origin='lower');
        axes[3].plot(self.r,mx_mean,self.r,my_mean)
        axes[3].legend(['<-$r \\Sigma v_r$>', '-<$\\Sigma v_\\phi$>'],loc='best')
        axes[2].plot(self.r,divm_mean,self.r,drmx_mean,self.r,dpmy_mean)
        axes[2].legend(['< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >', '-<$\\nabla_r (r \\Sigma v_r)$>', '-<$\\nabla_\\phi(\\Sigma v_\\phi)$>'],loc='best')


    def calculate_total_mdot(self,logspacing=True):
        pp,rr = meshgrid(self.phi,self.r)
        if logspacing:
            dr = diff(log(self.r))[0]
            norm = rr
        else:
            dr = diff(self.r)[0]
            norm = 1

        dp = diff(self.phi)[0]

        dens = self.dens.data
        vp = self.vpc
        vr = self.vr.data

        mx = -rr * dens*vr
        my = -dens*vp

        drmx,_ = gradient(mx,dlr,dp)
        _,dpmy = gradient(my,dlr,dp)

        drmx /= (rr*norm)
        dpmy /= rr
        divm = drmx + dpmy
        return dens,divm, divm.mean(axis=1),drmx.mean(axis=1), dpmy.mean(axis=1), mx.mean(axis=1), my.mean(axis=1)



    def summary(self):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)
    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('$<v_r$>',fontsize=20)
    	axd.set_ylabel('$<\\Sigma>$',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)
    #		axm.set_yscale('symlog',linthreshy=1e-7)
    #		axv.set_ylim(-.001,.001)
    	axd.set_title('t = %.1f = %.1f P = %.1f t_visc' % (self.t,self.t/(2*pi*self.a**(1.5)),self.t/self.tvisc))
    def streams(self,rlims=None,plims=None,ax=None,noise=.1,clrbar=True,**kargs):
        draw_flag = False
        if ax == None:
            fig=figure()
            ax=fig.add_subplot(111)
            draw_flag = True

        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)

        rinds = (self.r<=rlims[1])&(self.r>=rlims[0])
        pinds = (self.phi<=plims[1])&(self.phi>=plims[0])

        vr = self.vr.data[rinds,:][:,pinds]
        vp = self.vp.data[rinds,:][:,pinds]
        dens = self.dens.data[rinds,:][:,pinds]
        rr,pp = meshgrid(self.r[rinds],self.phi[pinds])
        lr = log(self.r[rinds])
        phi = self.phi[pinds]
        cmap = kargs.pop('cmap',viridis)

        line2d= ax.pcolormesh(log(rr),pp,log10(dens.transpose()),cmap=cmap,shading='gouraud')
        if clrbar:
            cbar = colorbar(line2d,ax=ax)
            cbar.set_label('$\\log_{10}{\\Sigma}$',fontsize=20)
        ax.streamplot(log(self.r[rinds]),self.phi[pinds],vr.transpose(),vp.transpose(),**kargs)
        ax.set_xlabel('$\ln(r)$',fontsize=20)
        ax.set_ylabel('$\phi$',fontsize=20)
        ax.set_ylim(plims)
        ax.set_xlim((log(rlims[0]),log(rlims[1])))
        ax.axvline(log(self.a+2*self.rh),color='k',linewidth=3)
        ax.axvline(log(self.a-2*self.rh),color='k',linewidth=3)


        rh,ph = self.calculate_circle(self.rh,rlims)
        rs,ps = self.calculate_circle(self.dens.aspectratio*self.a*self.dens.thicknesssmoothing,rlims)

        ax.plot(log(rh), ph,'-r',linewidth=3)
        ax.plot(log(rs),ps,'--r',linewidth=3)

        sep_lines = self.separatrix(lr,phi,vr.transpose(),vp.transpose(),noise=noise,npoints=10)
        for line in sep_lines:
            ax.plot(line[:,0],line[:,1],'-w',linewidth=2)
        if draw_flag:
            fig.canvas.draw()
#        return log(self.r[rinds]),self.phi[pinds],rr,pp,vr.transpose(),vp.transpose()

    def calculate_circle(self,d,rlims=None):
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)

        rinds = (self.r>=rlims[0])&(self.r<=rlims[1])
        newinds = (self.r[rinds]>=(self.a-d))&(self.r[rinds]<=(self.a+d))
        circle_upper = arccos((self.r[rinds][newinds]**2 + self.a**2 - d**2)/(2*self.a*self.r[rinds][newinds]))
        circle_lower=  -circle_upper

        circle_r = hstack((array([self.a-d]),self.r[rinds][newinds]))
        circle_r = hstack((circle_r,array([self.a+d])))
        circle_r = hstack((circle_r,self.r[rinds][newinds][::-1]))

        circle_phi = hstack((array([0]),circle_upper))
        circle_phi = hstack((circle_phi,array([0])))
        circle_phi = hstack((circle_phi,circle_lower[::-1]))
        return circle_r, circle_phi

    def stagnation(self,lr,phi,vr,vp):
        cr = contour(lr,phi,vr,levels=(0,),colors='w')
        cp = contour(lr,phi,vp,levels=(0,),colors='w')

        stag_points = []
        for pathr in cr.collections[0].get_paths():
            vertr = pathr.vertices
            segsr =[ vstack( (vertr[i,:],vertr[i+1,:])) for i in range(vertr.shape[0]-1)]
            for pathp in cp.collections[0].get_paths():
                vertp = pathp.vertices
                segsp =[ vstack( (vertp[i,:],vertp[i+1,:])) for i in range(vertp.shape[0]-1)]
                for sr in segsr:
                    for sp in segsp:
                        if intersect(sr,sp):
                            stag_points.append(get_intersection(sr,sp))

        return stag_points

    def separatrix(self,lr,phi,vr,vp,noise=0.1,npoints=10):
        ivp = interp2d(lr,phi,vp)
        ivr = interp2d(lr,phi,vr)
        dlr = lr[1]-lr[0]
        dp = phi[1]-phi[0]

        plt.figure()
        stag_points = self.stagnation(lr,phi,vr,vp)
        plt.close()
        uvals=np.array([x[0] for x in stag_points])
        inds = np.argsort(uvals)
#        for x in stag_points:
#            plot(x[0],x[1],'ow',markersize=20)
        sL = stag_points[inds[0]]
        sR = stag_points[inds[-1]]
        lines=[]

        tempx = noise*np.sqrt(dlr**2 + dp**2)

        for i in range(npoints*2):
            theta = 2*np.pi *  i/(float(2*npoints))
            u0 = sL[0] + tempx * np.cos(theta)
            p0 = sL[1] + tempx*np.sin(theta)

#            u0 = stag_points[0][0] + dlr*noise*2*(-.5 + rand())
#            p0 = stag_points[0][1] + dp*noise*2*(-.5+rand())
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,False))
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,True))
            u0 = sR[0] + tempx * np.cos(theta)
            p0 = sR[1] + tempx*np.sin(theta)

#            u0 = stag_points[0][0] + dlr*noise*2*(-.5 + rand())
#            p0 = stag_points[0][1] + dp*noise*2*(-.5+rand())
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,False))
            lines.append(self.get_stream(u0,p0,ivp,ivr,lr,phi,True))


        return lines
    def get_stream(self,u0,p0,ivp,ivr,lr,phi,reverse=False):
        dlr = lr[1]-lr[0]
        dp = phi[1]-phi[0]

        lr_min,lr_max = (lr.min(),lr.max())
        phi_min,phi_max = (phi.min(),phi.max())

        l = np.min((dlr,dp))
        nmax=1000

        res = np.array([u0,p0])
        uj = u0
        pj=  p0
        breakflag = False
        for j in range(nmax):
            dx, dy = self.euler_step(uj,pj,l,ivp,ivr,reverse)
            uj += dx
            pj += dy
            if uj >= lr_max:
                uj = lr_max
#                print 'Stream lnr > lnr_max'
                breakflag = True
            if uj <= lr_min:
 #               print 'Stream lnr < lnr_min'
                uj = lr_min
                breakflag = True
            if pj >= phi_max:
  #              print 'Stream phi > phi_max'
                pj = phi_max
                breakflag = True
            if pj <= phi_min:
   #             print 'Stream phi < phi_min'
                pj = phi_min
                breakflag = True
            if breakflag:
                break
            else:
               res=np.vstack( (res,np.array([uj,pj])))
    #    print 'Iterated past max iterations of %d' % nmax
        return res

    def euler_step(self,u0,p0,l,ivp,ivr,reverse=False):
        safety_fac = self.safety_fac
        h = 1.0
        if reverse:
            h = -1.0

        vp_val = ivp(u0,p0)[0]
        vr_val = ivr(u0,p0)[0]

        h *= safety_fac*l/np.sqrt( vp_val**2 + vr_val**2)
        return h*vr_val, h*vp_val



    def animate(self,irange):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)

        s0 = Sim(0,directory=self.directory)
        axd.plot(s0.r,s0.dbar,'--k',linewidth=3)
    	axv.plot(s0.r,s0.vrbar,'--k',linewidth=3)
    	axm.plot(s0.r,s0.mdot,'--k',linewidth=3)
        axv.set_ylim(-.05,.05)

    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('<$v_r$>',fontsize=20)
    	axd.set_ylabel('<$\\Sigma$>',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)

    	dbar = zeros((len(self.dbar),len(irange)))
    	vrbar = zeros((len(self.dbar),len(irange)))
    	mdot = zeros(dbar.shape)
    	t = zeros(len(irange))
    	a = zeros(len(irange))

    	t[0] = self.t
    	a[0] = self.a
    	dbar[:,0] = self.dbar
    	vrbar[:,0] = self.vrbar
    	mdot[:,0] = self.mdot
    	for i,j in enumerate(irange[1:],start=1):
            self.__init__(j,directory=self.directory)
            t[i] = self.t
            dbar[:,i] = self.dbar
            vrbar[:,i] = self.vrbar
            mdot[:,i] = self.mdot
            a[i] = self.a

#    	axv.set_ylim(-.001,.001)
#    	axm.set_ylim((mdot.min(),mdot.max()))
#    	axm.set_yscale('symlog',linthreshy=1e-8)
    	for i in range(len(t)):
    		lined.set_ydata(dbar[:,i])
    		linev.set_ydata(vrbar[:,i])
    		linem.set_ydata(mdot[:,i])
    		axd.set_title('t = %.1f = %.1f P' % (t[i],t[i]/(2*pi*a[i]**(1.5))),color='w')
    		fig.canvas.draw()


    def steady_state_calc(self):
        self.mdot0_ss = 3*pi*self.nu0*(self.dens.ymax*self.dens.sigmaout-self.dens.ymin*self.dens.sigmain)/(np.sqrt(self.dens.ymax)-np.sqrt(self.dens.ymin))
        self.lam0_ss = 2*pi*self.dens.sigmain*self.dens.ymax + 2*self.mdot0_ss*(np.sqrt(self.r)-np.sqrt(self.dens.ymin))/(3*self.nu0)
        self.vr0_ss = -self.mdot0_ss/self.lam0_ss
        self.dbar0_ss = self.lam0_ss/(2*pi*self.r)


    def refine_grid(self,levels=3,logspacing=True,fname='domain_y.dat',save=False,savedir='./'):
        if levels > 3:
            print 'Too many levels, maximum is 3'
            return
        rad = loadtxt(fname)
        hs_zone = lambda r: abs(r-self.a) <= self.rh*2
        rh_zone = lambda r: abs(r-self.a) <= self.rh
        eps = self.dens.thicknesssmoothing * self.dens.aspectratio * self.a
        soft_zone = lambda r: abs(r-self.a) <= eps

        print 'Planet at %f\nHs zone at %f\nHill zone at %f\nSoft zone at %f' %(self.a,2*self.rh,self.rh,eps)

        if eps > self.rh:
            if eps > 2*self.rh:
                ind_func = [soft_zone,hs_zone,rh_zone]
            else:
                ind_func = [hs_zone,soft_zone,rh_zone]
            print 'smoothing zone larger than hill zone'
        else:
            ind_func = [hs_zone,rh_zone,soft_zone]

        for j in range(levels):
            ind = ind_func[j](rad)
            ri = rad[ind][0]
            ro = rad[ind][-1]

            if logspacing:
                spacing = self.dlr * 2.**(-(j+1))
                rad1 = np.exp(np.linspace(np.log(ri),np.log(ro),np.log(ro/ri)/spacing+1))
            else:
                spacing = self.dr[0] * 2.**(-(j+1))
                rad1 = np.linspace(ri,ro,(ro-ri)/spacing+1)

            leftr = rad[ rad < ri]
            rightr = rad[ rad > ro]

            rad = np.hstack( (leftr,rad1,rightr) )

        print 'Went from %d points to %d points' % (self.dens.ny,len(rad)-7)

        if save:
            call(['cp',fname,fname+'.backup'])
            with open(savedir+fname,'w') as f:
                f.write('\n'.join(rad.astype(str)))

        return rad
    def refine_phi(self,j=2,fname='domain_x.dat',save=False,savedir='./'):
        x = np.loadtxt(fname)
        dx = diff(x)[0]
        x = np.linspace(x[0],x[-1],(x[-1]-x[0])/(dx/j)+1)


        print 'Went from %d points to %d points' % (self.dens.ny,len(x)-1)

        if save:
            call(['cp',fname,fname+'.backup'])
            with open(savedir+fname,'w') as f:
                f.write('\n'.join(x.astype(str)))

        return x
    def linear_torque(self,r,eps=.2,c=2./3):
        norm = eps*pi*self.mp**2
        ha = self.a*self.dens.aspectratio
        x = (r-self.a)/ha

        sgn = (x>=0).astype(int)
        sgn[ sgn == 0 ] = -1

        indR = x >=c
        indL = x <= -c
        res = zeros(r.shape)
        res[indR] = (self.a/(r[indR]-self.a))**4
        res[indL] = (r[indL]/(r[indL]-self.a))**4

        return res * norm * sgn


    def load_torques(self):
        try:
            dat=fromfile('torq_1d_Y_raw_planet_0.dat')
            nt = len(dat)/self.dens.ny
        except IOError:
            print "Can't find torque profile file"
            dat = np.zeros((1,self.dens.ny))
            nt = 1

        try:
            tq=loadtxt('torq_planet_0.dat')
        except IOError:
            print "Can't find total torque file"
            tq = np.zeros((nt,2))

        try:
            mass=fromfile('mass_1d_Y_raw.dat')
            mass = mass.reshape(nt,self.dens.ny)
        except IOError:
            print "Can't find mass profile"
            mass = np.ones((nt,self.dens.ny))


        dat = dat.reshape(nt,self.dens.ny)
        tq[:,-1] *= -self.mp
        dat *= -self.mp


        self.int_lambda = zeros(dat.shape)
        self.int_lambda[:,1:] = cumtrapz(dat,axis=1)
        for i in range(nt):
            dat[i,:] *= self.dphi[0]*self.dens.ymed/self.vol[:,0]
            mass[i,:] *= self.dphi[0]/self.vol[:,0]
            mass[i,:] /= 2*pi
        self.total_torque = tq
        self.lambda_prof = dat
        self.mass_prof = mass
        self.torque_dens_prof = dat/mass
        return
    def compute_fft(self):
        import numpy.fft as ft
        dhat = ft.rfft(self.dens.data,axis=1)/(self.dens.nx-1)
        vrhat = ft.rfft(self.vr.data,axis=1)/(self.dens.nx-1)
        vphat = ft.rfft(self.vp.data,axis=1)/(self.dens.nx-1)
        p_dhat = trapz(dhat*conj(dhat)*self.dr[:,np.newaxis],axis=0)
        p_vrhat = trapz(vrhat*conj(vrhat)*self.dr[:,np.newaxis],axis=0)
        p_vphat = trapz(vphat*conj(vphat)*self.dr[:,np.newaxis],axis=0)
        return dhat,vrhat,vphat,p_dhat,p_vrhat,p_vphat

def ccw(a,b,c):
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])

def intersect(segment_1, segment_2):
    a1 = segment_1[0,:]
    a2 = segment_1[1,:]
    b1 = segment_2[0,:]
    b2 = segment_2[1,:]
    return ccw(a1,b1,b2) != ccw(a2,b1,b2) and ccw(a1,a2,b1) != ccw(a1,a2,b2)
def get_intersection(segment_1, segment_2):
    if intersect(segment_1,segment_2):
#        print 'Calculating Intersection'
        m1 = segment_1[1,:]-segment_1[0,:]
        m1 = m1[1]/m1[0]
        m2 = segment_2[1,:] - segment_2[0,:]
        m2 = m2[1]/m2[0]
        lhs = np.array([[-m1,1.],[-m2,1.]])
        rhs = np.array([[ segment_1[0,1] - m1*segment_1[0,0]],[segment_2[0,1]-m2*segment_2[0,0]]])
        return dot( np.linalg.inv(lhs),rhs).reshape(2)
    else:
#        print 'Lines do not intersect'
        return False
def boundary_progress(i,fig=None,axes=None,j=0,alpha=0.01,h=0.05,flaringindex=0):
    config_flag = True
    if fig is None:
        fig,axes=subplots(3,1,sharex=True)
        config_flag = False
    try:
        rho = fromfile('gasdens{0:d}_0.dat'.format(i))
        vr = fromfile('gasvy{0:d}_0.dat'.format(i))
        ym = loadtxt('domain_y.dat')
        xm = loadtxt('domain_x.dat')
        ymed = (ym[1:] + ym[:-1])/2
        xmed = (xm[1:] + xm[:-1])/2
        ny = len(ymed)-6
        nx = len(xmed)
        vr = vr.reshape(ny+6,nx)
        rho = rho.reshape(ny+6,nx)
    except IOError:
        rho = fromfile('gasdens{0:d}.dat'.format(i))
        vr = fromfile('gasvy{0:d}.dat'.format(i))
        ym = loadtxt('domain_y.dat')
        xm = loadtxt('domain_x.dat')
        ymed = (ym[1:] + ym[:-1])/2
        xmed = (xm[1:] + xm[:-1])/2
        ny = len(ymed)-6
        nx = len(xmed)
        rho = rho.reshape(ny,nx)
        vr = vr.reshape(ny,nx)
#    mdot = 2*pi*ymed * (rho*vr).mean(axis=1)
    #mdot = 2*ymed[:-1]*pi*(rho[:-1,:]*((vr[:1,:]+vr[:-1,:])/2)).mean(axis=1)
    rho = rho.mean(axis=1)
    vr = vr.mean(axis=1)
    nu = alpha*h*h*ymed**(2*flaringindex+0.5)
    fj = 3*pi*nu*rho*sqrt(ymed)

    axes[0].plot(ymed,rho,'.-')
    axes[1].plot(ymed,fj,'.-')
    axes[2].plot(ym[:-1],vr,'.-')

    axes[0].set_ylabel('$\\Sigma$',fontsize=20)
    axes[1].set_ylabel('$F_\\nu$',fontsize=20)
    axes[2].set_ylabel('$v_r$',fontsize=20)


    axes[2].set_xlabel('$r$',fontsize=20)
#    axes[2].plot(ymed,mdot)

    if not config_flag:
        for ax in axes:
            ax.minorticks_on()
#            ax.set_xscale('log')
#        axes[0].set_yscale('log')
        subplots_adjust(hspace=0)
    return fig,axes,ymed,ym,rho,fj,vr



def time_avg_mdot(irange):

    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])
    ymed = ymed[:-1]
    fac = 2*pi*ymed
    facm = 2*pi*ym[1:-1]
    res = zeros((len(ymed),len(irange)))
    res_rho = zeros((len(ymed),len(irange)))
    for i,t in enumerate(irange):
        rho,momy,_,_,_,momys = load_single_time_data(t)
        res[:,i] = (momys).mean(axis=1)
        res_rho[:,i] = rho.mean(axis=1)

    return -facm*res.mean(axis=1), -fac*res.std(axis=1),res_rho.mean(axis=1)*fac,res_rho.std(axis=1)*fac,ymed



def load_single_time_data(i):
    vx,vy,rho = load_single_time(i)

    vxc = 0.5*(vx.data[:,:-1] + vx.data[:,1:])
    vyc = 0.5*(vy.data[:-1,:] + vy.data[1:,:])

    momy = vyc*rho.data[:-1,:]
    momx = vxc*rho.data[:,:-1]

    momys = .5*(rho.data[1:,:] + rho.data[:-1,:])*vy.data[1:,:]



    return rho.data[:-1,:],momy,vxc,vyc,rho.ymed[:-1],momys

def vortencity(rho,vx,vy):
    Tx,Rx = meshgrid(vx.x,vx.y)
    Ty,Ry = meshgrid(vy.x,vy.y)
    rvx = Rx*(vx.data)

    curldata = (( rvx[1:,1:]-rvx[:-1,1:])/(Rx[1:,1:]-Rx[:-1,1:]) -
            (vy.data[1:,1:] - vy.data[1:,:-1])/(Ty[1:,1:]-Ty[1:,:-1]))

    curl = copy.deepcopy(vx)
    curl.nx = curl.nx-1
    curl.ny = curl.ny-1
    curl.x = vx.x[:-1]
    curl.y = vy.y[:-1]

    rho_corner = .25*(rho.data[1:,1:] + rho.data[:-1,:-1] + rho.data[1:,:-1] + rho.data[:-1,1:])

    T,R = meshgrid(curl.x,curl.y)
    curl.data = (curldata/R + 2*rho.omegaframe)/rho_corner
    return curl

def load_single_time(i):

    rho = Field('gasdens{0:d}.dat'.format(i))
    vx = Field('gasvx{0:d}.dat'.format(i),staggered='x')
    vy = Field('gasvy{0:d}.dat'.format(i),staggered='y')


    return vx,vy,rho

def load_flux(i):

    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])

    xm = loadtxt('domain_x.dat')
    xmed = .5*(xm[1:] + xm[:-1])
    dx = diff(xm)
    dy = diff(ym)

    ny = len(ymed)
    nx  = len(xmed)


    rho = fromfile('gasdens{0:d}.dat'.format(i)).reshape(ny,nx)
    vy = fromfile('gasvy{0:d}.dat'.format(i)).reshape(ny,nx)

    dqm = zeros(rho.shape)
    dqp = zeros(rho.shape)
    slope = zeros(rho.shape)
    mdot = zeros(rho.shape)

    surfy = outer(ym,dx)
    dxx,dyy = meshgrid(dx,dy)
    dqm = zeros(rho.shape)
    dqp = zeros(rho.shape)
    slope1 = zeros(rho.shape)
    mdot1=zeros(rho.shape)
    surfy1 = outer(ym,dx)
    dqm[1:-1] = (rho[1:-1,:]-rho[:-2,:])/dyy[1:-1,:]
    dqp[1:-1] = (rho[2:,:]-rho[1:-1,:])/dyy[2:,:]
    slope1[1:-1,:] = ( (dqm*dqp)[1:-1,:] > 0).astype(int) * (2*(dqp*dqm)/(dqp+dqm))[1:-1,:]
    mdot[1:-1,:] = (vy[1:-1,:]>0).astype(int) * vy[1:-1,:]*(rho[:-2,:]+.5*slope[:-2,:]*dyy[:-2,:])*surfy1[1:-2,:]
    mdot[1:-1,:] += (vy[1:-1,:]<=0).astype(int) * vy[1:-1,:]*(rho[1:-1,:]-.5*slope[1:-1,:]*dyy[1:-1,:])*surfy1[1:-2,:]

    return ym[1:-2],-mdot.sum(axis=1)[1:-1]

def calc_mass_flux(i,cfl=0.5,h=0.05,dt=None,fac=1,plot_flag=False):


    ym = loadtxt('domain_y.dat')[3:-3]
    ymed = .5*(ym[1:] + ym[:-1])

    xm = loadtxt('domain_x.dat')
    xmed = .5*(xm[1:] + xm[:-1])
    dx = diff(xm)
    dy = diff(ym)

    ny = len(ymed)
    nx  = len(xmed)

    rho = fromfile('gasdens{0:d}.dat'.format(i)).reshape(ny,nx)
    vx = fromfile('gasvx{0:d}.dat'.format(i)).reshape(ny,nx)
    vy = fromfile('gasvy{0:d}.dat'.format(i)).reshape(ny,nx)
    rho1 = fromfile('gasdens{0:d}.dat'.format(i+1)).reshape(ny,nx)
    vy1 = fromfile('gasvy{0:d}.dat'.format(i+1)).reshape(ny,nx)



    flux_simple = zeros(rho.shape)
    rhonew_simple = copy.copy(rho)



    if dt is None:
        dt = cfl * min( dy/(h/sqrt(ymed)))

    print 'Using dt = %f' % dt

    rhostar = zeros(rho.shape)
    slope = zeros(rho.shape)

    zone_size_y = lambda j: ym[j+1]-ym[j]

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            dqm = (rho[i,j]-rho[i-1,j]) / zone_size_y(i)
            dqp = (rho[i+1,j]-rho[i,j]) / zone_size_y(i+1)
            if dqm*dqp <= 0:
                slope[i,j] = 0
            else:
                slope[i,j] = 2*dqp*dqm/(dqp+dqm)


 #   dqm = zeros(rho.shape)
 #   dqp = zeros(rho.shape)
 #   slope = zeros(rho.shape)
 #   mdot = zeros(rho.shape)

 #   surfy = outer(dy,dx)
 #   dxx,dyy = meshgrid(dx,dy)
 #   vol = outer(.5*(ym[1:]**2-ym[:-1]**2),dx)

 #   dqm[1:-1,:] =(rho[1:-1,:]  - rho[:-2,:])/dy[1:-1,:]
 #   dqp[1:-1,:] =(rho[2:,:]  - rho[1:-1,:])/dy[2:,:]
 #   ind = sign(dqm*dqp)
 #   slope = (ind>0).astype(int) * 2*dqm*dqp/(dqm+dqp)
 #   mdot[1:-1,:] = (vy>0).astype(int)*(rho[:-2,:]+.5*slope[:-2,:]*dy[:-2,:]) + (vy<=0).astype(int)*(rho[1:,:]-.5*slope[1:,:]*dy[1:,:])
 #   mdot[1:-1,:] *= surfy[1:-1,:]*vy[1:-1,:]



    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            if vy[i,j] > 0:
                rhostar[i,j] = rho[i-1,j] + .5*slope[i-1,j]*( zone_size_y(i-1) - vy[i,j]*dt)
            else:
                rhostar[i,j] = rho[i,j] - .5*slope[i,j]*( zone_size_y(i) + vy[i,j]*dt)


    rhonew = copy.copy(rho)
    flux = zeros(rho.shape)
    vol = (2*pi/ny)*.5*(ym[1:]**2 - ym[:-1]**2)
    surfy = (2*pi/ny)*ym
    inv = 1./vol

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            flux[i,j] = -(vy[i,j]*rhostar[i,j]*surfy[i]-vy[i+1,j]*rhostar[i+1,j]*surfy[i+1])
            rhonew[i,j] -= dt*inv[i]*flux[i,j]

    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            flux_simple[i,j] = surfy[i+1]*.5*(rho[i,j] + rho[i+1,j])*vy[i+1,j] - surfy[i]*.5*(rho[i-1,j]+rho[i,j])*vy[i,j]
            rhonew_simple[i,j] += - dt*inv[i] * flux_simple[i,j]


    fluxSp = zeros(rho.shape)
    fluxSm = zeros(rho.shape)
    fluxp = zeros(rho.shape)
    fluxm = zeros(rho.shape)
    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            fluxSp[i,j] = surfy[i+1]*.5*(rho[i,j] + rho[i+1,j])*vy[i+1,j]
            fluxSm[i,j] = surfy[i]*.5*(rho[i-1,j]+rho[i,j])*vy[i,j]
            fluxp[i,j] = vy[i+1,j]*rhostar[i+1,j]*surfy[i+1]
            fluxm[i,j] = vy[i,j]*rhostar[i,j]*surfy[i]


    indt = (vy>0).astype(int)
    indf = (~(vy>0)).astype(int)

    pi_p = zeros(rho.shape)
    pi_m = zeros(rho.shape)
    fluxup = zeros(rho.shape)
    for i,yc in enumerate(ymed[1:-1],start=1):
        for j,xc in enumerate(xmed):
            pi_p[i,j] = vy[i,j] * rho[i,j]*surfy[i]
            pi_m[i,j] = vy[i,j] * rho[i-1,j]*surfy[i]
            ap = slope[i-1,j]
            am = slope[i,j]
            if vy[i,j] > 0:
                fluxup[i,j] = vy[i,j]*(rho[i-1,j] + .5*ap*zone_size_y(i-1) ) * surfy[i]
            else:
                fluxup[i,j] = vy[i,j]* (rho[i,j] - .5*am*zone_size_y(i) )*surfy[i]

    if plot_flag:
        fig,axes=subplots(5,1,sharex=True)
        axes[0].set_title('$\\Delta t = %f$'%dt,fontsize=20)


        axes[4].plot(ym[:-1],vy.mean(axis=1),'-b')
        axes[4].set_ylabel('$v_r$',fontsize=20)

        axes[3].plot(ym[1:][1:-1],fluxp.mean(axis=1)[1:-1],'-b')
        axes[3].plot(ym[:-1][1:-1],fluxm.mean(axis=1)[1:-1],'--b')
        axes[3].plot(ym[1:][1:-1],fluxSp.mean(axis=1)[1:-1],'-r')
        axes[3].plot(ym[:-1][1:-1],fluxSm.mean(axis=1)[1:-1],'--r')
    #    axes[3].plot(ymed[:-1][1:-1],(ymed[:-1]*(2*pi/nx)* (.5*(vy[1:,:]+vy[:-1,:])*rho[:-1,:]).mean(axis=1))[1:-1],'--r')
    #    axes[3].plot(ym[1:-2],pi_p.mean(axis=1)[1:-1],'-k')
    #    axes[3].plot(ym[1:-2],pi_m.mean(axis=1)[1:-1],'--k')
        axes[3].plot(ym[1:-2],fluxup.mean(axis=1)[1:-1],'-m')
        axes[3].set_ylabel('$F$',fontsize=20)

        axes[2].plot(ymed[1:-1], 1.e-16 +abs( ( (rhonew-rho1)/rho1 ).mean(axis=1))[1:-1],'-b')
        axes[2].plot(ymed[1:-1], 1.e-16 +abs( ( (rhonew_simple-rho1)/rho1 ).mean(axis=1))[1:-1],'-r')
        axes[2].set_yscale('log')
        axes[2].set_ylabel('Relative Error',fontsize=15)

        axes[0].plot(ymed[1:-1],(rhonew-rho).mean(axis=1)[1:-1],'-b')
        axes[0].plot(ymed[1:-1],(rhonew_simple-rho).mean(axis=1)[1:-1],'-r')
        axes[0].plot(ymed[1:-1],(rho1-rho).mean(axis=1)[1:-1]*fac,'-g')
        axes[0].set_ylabel('$\\Delta\\Sigma$',fontsize=20)
        axes[0].legend(['Mock step', '$\\Sigma^*=\\Sigma_{avg}$','FARGO'],loc='upper right')


        axes[1].plot(ymed[1:-1],flux.mean(axis=1)[1:-1],'-b')
        axes[1].plot(ymed[1:-1],flux_simple.mean(axis=1)[1:-1],'-r')
        axes[1].set_ylabel('$\\Delta F$',fontsize=20)

        axes[-1].set_xlabel('$r$',fontsize=20)
        for ax in axes:
            ax.minorticks_on()

    return ym[:-1][1:-1],-fluxup.sum(axis=1)[1:-1]

def check_mdot():
    t = loadtxt('mass.dat')
    rho = fromfile('mass_1d_Y_raw.dat')
    momy = fromfile('momy_1d_Y_raw.dat')

    ym = loadtxt('domain_y.dat')[3:-3]
    xm = loadtxt('domain_x.dat')

    ymed = .5*(ym[:-1]+ym[1:])
    xmed = .5*(xm[:-1] + xm[1:])

    surf = .5*(ym[1:]**2 - ym[:-1]**2)
    vol = surf * diff(xm)[0] # Uniform gird in phi

    dy = diff(ym)
    dly = dy/ymed

    ny = len(ymed)
    nx = len(xmed)
    nt = len(rho)/ny

    if len(rho)/ny != len(momy)/ny or len(rho)/ny != len(t[:,0]) or len(momy)/ny != len(t[:,0]):
        print "File's loaded at different times! Try again."
        return

    rho = rho.reshape(nt,ny)
    momy = momy.reshape(nt,ny)
    t = t[:,0]


    for i in range(rho.shape[0]):
        rho[i,:] /= vol
        momy[i,:] /= vol
        momy[i,:] *= -2*pi*ymed
    rho /= nx

    momy /= nx



    return ymed,rho,momy,t


def check_mdot_single(i):
    rho = fromfile('gasdens{0:d}.dat'.format(i))
    vy = fromfile('gasvy{0:d}.dat'.format(i))

    rho0 = fromfile('gasdens{0:d}.dat'.format(i-1))

    ym = loadtxt('domain_y.dat')[3:-3]
    xm = loadtxt('domain_x.dat')

    ymed = .5*(ym[:-1]+ym[1:])
    xmed = .5*(xm[:-1] + xm[1:])

    surf = .5*(ym[1:]**2 - ym[:-1]**2)
    vol = surf * diff(xm)[0] # Uniform gird in phi

    dy = diff(ym)
    dly = dy/ymed

    ny = len(ymed)
    nx = len(xmed)

    rho = rho.reshape(ny,nx)
    rho0 = rho0.reshape(ny,nx)
    vy = vy.reshape(ny,nx)


 #   dbar = rho.mean(axis=1)

#    momy = -2*pi*ymed[:-1]*(rho[:-1,:]*.5*(vy[1:,:] + vy[:-1,:])).mean(axis=1)


#    momys = -2*pi*ym[1:-1]*(  .5*(rho[1:,:] + rho[:-1,:])*vy[1:,:] ).mean(axis=1)


 #   rho_1 = interp1d(ymed,rho,axis=0,kind='slinear')
 #   rho_2 = interp1d(ymed,rho,axis=0,kind='quadratic')
    rho_3 = interp1d(ymed,rho,axis=0,kind='cubic')


#    momy1 =-2*pi*ym[1:-1]* (rho_1(ym[1:-1])*vy[1:,:]).mean(axis=1)
#    momy2 =-2*pi*ym[1:-1]* (rho_2(ym[1:-1])*vy[1:,:]).mean(axis=1)
    momy3 =-2*pi*ym[1:-1]* (rho_3(ym[1:-1])*vy[1:,:]).mean(axis=1)


    params = Parameters()

    dt = params.dt * params.ninterm

    dlamdt = 2*pi*ymed*(rho-rho0).mean(axis=1)/dt


    figure()
    plot(ymed,dlamdt*dy)
    figure();
    plot(ym[1:-2],diff(momy3))



    return ym[1:-1],momy3



def grab_axis(i):
    return figure(i).get_axes()[0]


