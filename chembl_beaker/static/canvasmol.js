/*****************************************************************************************/
// Globals
/*****************************************************************************************/
var INIT_NAMES = ["c60.sdf"];
//var INIT_NAMES = ["axis.mol"];

var MENU_MAP = {
"Simple"    :0,
"Anandamide":"anandamide.sdf",
"Aspirin"   :"aspirin.sdf",
"ATP"       :"atp.sdf",
"Caffeine"  :"caffeine.sdf",
"Cantharidin":"cantharidin.sdf",
"Capsaicin" :"capsaicin.sdf",
"Carotene"  :"carotene.mol",
"Chlorophyll":"chlorophyll.mol",
"Cholesterol":"cholesterol.sdf",
"Cocaine"   :"cocaine.sdf",
"Cubane"    :"cubane.sdf",
"Curcumin"  :"curcumin.sdf",
"DMT"       :"dmt.sdf",
"Ethanol"   :"ethanol.pdb",
"Fructose"  :"fructose.sdf",
"Glucose"   :"glucose.sdf",
"LSD"       :"lsd.sdf",
"Lycopene"  :"lycopene.mol",
"MDMA"      :"mdma.sdf",
"Modafinil" :"modafinil.sdf",
"Morphine"  :"morphine.sdf",
"Nicotine"  :"nicotine.sdf",
"Oxytocin"  :"oxytocin.sdf",
"PCP"       :"pcp.sdf",
"Penicillin":"penicillin.sdf",
"Psilocybin":"psilocybin.sdf",
"Sildenafil":"sildenafil.sdf",
"Strychnine":"strychnine.sdf",
"Teflon"    :"ptfe.sdf",
"THC"       :"thc.sdf",
"TNT"       :"tnt.sdf",

"Crystals"  :0,
"Aluminium oxide":"Al2O3.pdb",
"Copper"    :"cu.pdb",
"Fluorite"  :"caf2.pdb",
"Salt"      :"nacl.pdb",
"YBCO superconductor":"ybco.pdb",

"Carbon"    :0,
"Diamond"   :"diamond.pdb",
"Graphene"  :"graphene.sdf",
"Graphite"  :"graphite.pdb",
"Buckyball" :"c60.sdf",
"Buckytube" :"nanotube.sdf",

"Bio"       :0,
"Daptomycin"        :"daptomycin.pdb",
"DNA"               :"dna.sdf",
"DNA crystal"       :"3GBI.pdb",
"Transfer RNA"      :"4TNA.pdb",
"Adenovirus"        :"adenovirus.pdb",
"Black beetle virus"  :"2bbv.pdb",
"Tobacco mosaic virus":"2OM3.pdb",
"Cyclosporin complex":"2X2C.pdb",
"KaiC protein"       :"2GBL.pdb",
"Human prion protein":"1OEI.pdb",
"Cadherin"           :"1L3W.pdb",
"Taq polymerase"     :"2GHO.pdb"
};

var MENU_MAP_LO = {};
for(var i in MENU_MAP)
    if(MENU_MAP[i]) MENU_MAP_LO[i.toLowerCase()] = MENU_MAP[i];

var PARS = {
"default": {'draw_atoms'  :1,
         'draw_bonds'     :1,
         'auto_x'         :1,
         'auto_y'         :1,
         'auto_z'         :0,
         'bonds_autowidth':0,
         'bonds_gradient' :0,
         'atom_colors'    :1,
         'glow'           :1,
         'flat_shading'   :0,
         'sphere_shading' :1,
         'atom_radius'    :10
        },
"small": {'draw_atoms'    :1,
         'draw_bonds'     :1,
         'auto_x'         :1,
         'auto_y'         :1,
         'auto_z'         :0,
         'bonds_autowidth':0,
         'bonds_gradient' :0,
         'atom_colors'    :1,
         'glow'           :1,
         'flat_shading'   :0,
         'sphere_shading' :1,
         'atom_radius'    :5
        },
"full": {'draw_atoms'     :1,
         'draw_bonds'     :1,
         'auto_x'         :1,
         'auto_y'         :1,
         'auto_z'         :0,
         'bonds_autowidth':1,
         'bonds_gradient' :1,
         'atom_colors'    :1,
         'glow'           :1,
         'flat_shading'   :0,
         'sphere_shading' :1,
         'atom_radius'    :5
        },
"atoms":{'draw_atoms'     :1,
         'draw_bonds'     :0,
         'auto_x'         :1,
         'auto_y'         :1,
         'auto_z'         :0,
         'bonds_autowidth':0,
         'bonds_gradient' :0,
         'atom_colors'    :1,
         'glow'           :1,
         'flat_shading'   :0,
         'sphere_shading' :1,
         'atom_radius'    :5
        },
"atoms2":{'draw_atoms'     :1,
         'draw_bonds'     :0,
         'auto_x'         :0,
         'auto_y'         :0,
         'auto_z'         :0,
         'bonds_autowidth':0,
         'bonds_gradient' :0,
         'atom_colors'    :0,
         'glow'           :0,
         'flat_shading'   :1,
         'sphere_shading' :0,
         'atom_radius'    :5
        },
"bonds":{'draw_atoms'     :0,
         'draw_bonds'     :1,
         'auto_x'         :0,
         'auto_y'         :1,
         'auto_z'         :1,
         'bonds_autowidth':0,
         'bonds_gradient' :0,
         'atom_colors'    :0,
         'glow'           :0,
         'flat_shading'   :0,
         'sphere_shading' :1,
         'atom_radius'    :10
        }
};

var PAR_MAP = {
"2bbv"      :"atoms2",
"adenovirus":"atoms2",
"2GBL"      :"atoms2",
"2OM3"      :"atoms2",
"2X2C"      :"atoms2",
"2GHO"      :"atoms2",
"1OEI"      :"atoms2",
"1L3W"      :"atoms2",
"4TNA"      :"atoms",
"dna"       :"small",
"daptomycin":"small",
'carotene'  :"small",
'chlorophyll':"small",
'cholesterol':"small",
'oxytocin'  :"small",
"3GBI"      :"atoms",
"c60"       :"bonds",
"graphite"  :"bonds",
"diamond"   :"bonds",
"nanotube"  :"bonds"
};

function get_type(name) {
    var sname = name.substr(0, name.length-4);
    if(PAR_MAP[sname]==undefined)
        return "default";
    return PAR_MAP[sname];
}

// default colors
var BONDS_COLOR = [255,255,255];
var ATOMS_COLOR = [255,0,0];

var BACKGROUND_COLOR = [0,0,0];
var BACKGROUND_COLOR_RGB ='rgb('+BACKGROUND_COLOR[0]+','+BACKGROUND_COLOR[1]+','+BACKGROUND_COLOR[2]+')';

var CPK = {"h":[255,255,255],"he":[217,255,255],"li":[204,128,255],"be":[194,255,0],"b":[255,181,181],"c":[144,144,144],"n":[48,80,248],"o":[255,13,13],"f":[144,224,80],"ne":[179,227,245],"na":[171,92,242],"mg":[138,255,0],"al":[191,166,166],"si":[240,200,160],"p":[255,128,0],"s":[255,255,48],"cl":[31,240,31],"ar":[128,209,227],"k":[143,64,212],"ca":[61,255,0],"sc":[230,230,230],"ti":[191,194,199],"v":[166,166,171],"cr":[138,153,199],"mn":[156,122,199],"fe":[224,102,51],"co":[240,144,160],"ni":[80,208,80],"cu":[200,128,51],"zn":[125,128,176],"ga":[194,143,143],"ge":[102,143,143],"as":[189,128,227],"se":[255,161,0],"br":[166,41,41],"kr":[92,184,209],"rb":[112,46,176],"sr":[0,255,0],"y":[148,255,255],"zr":[148,224,224],"nb":[115,194,201],"mo":[84,181,181],"tc":[59,158,158],"ru":[36,143,143],"rh":[10,125,140],"pd":[0,105,133],"ag":[192,192,192],"cd":[255,217,143],"in":[166,117,115],"sn":[102,128,128],"sb":[158,99,181],"te":[212,122,0],"i":[148,0,148],"xe":[66,158,176],"cs":[87,23,143],"ba":[0,201,0],"la":[112,212,255],"ce":[255,255,199],"pr":[217,255,199],"nd":[199,255,199],"pm":[163,255,199],"sm":[143,255,199],"eu":[97,255,199],"gd":[69,255,199],"tb":[48,255,199],"dy":[31,255,199],"ho":[0,255,156],"er":[0,230,117],"tm":[0,212,82],"yb":[0,191,56],"lu":[0,171,36],"hf":[77,194,255],"ta":[77,166,255],"w":[33,148,214],"re":[38,125,171],"os":[38,102,150],"ir":[23,84,135],"pt":[208,208,224],"au":[255,209,35],"hg":[184,184,208],"tl":[166,84,77],"pb":[87,89,97],"bi":[158,79,181],"po":[171,92,0],"at":[117,79,69],"rn":[66,130,150],"fr":[66,0,102],"ra":[0,125,0],"ac":[112,171,250],"th":[0,186,255],"pa":[0,161,255],"u":[0,143,255],"np":[0,128,255],"pu":[0,107,255],"am":[84,92,242],"cm":[120,92,227],"bk":[138,79,227],"cf":[161,54,212],"es":[179,31,212],"fm":[179,31,186],"md":[179,13,166],"no":[189,13,135],"lr":[199,0,102],"rf":[204,0,89],"db":[209,0,79],"sg":[217,0,69],"bh":[224,0,56],"hs":[230,0,46],"mt":[235,0,38],
           "ds":[235,0,38],"rg":[235,0,38],"cn":[235,0,38],"uut":[235,0,38],"uuq":[235,0,38],"uup":[235,0,38],"uuh":[235,0,38],"uus":[235,0,38],"uuo":[235,0,38]};
var MENDELEEV = [
    //   1    2    3    4    5    6    7    8    9   10   11   12    13    14    15    16    17    18
    [  "h",   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,    0, "he" ],
    [ "li","be",   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  "b",  "c",  "n",  "o",  "f", "ne" ],
    [ "na","mg",   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, "al", "si",  "p",  "s", "cl", "ar" ],
    [  "k","ca","sc","ti", "v","cr","mn","fe","co","ni","cu","zn", "ga", "ge", "as", "se", "br", "kr" ],
    [ "rb","sr", "y","zr","nb","mo","tc","ru","rh","pd","ag","cd", "in", "sn", "sb", "te",  "i", "xe" ],
    [ "cs","ba","la","hf","ta", "w","re","os","ir","pt","au","hg", "tl", "pb", "bi", "po", "at", "rn" ],
    [ "fr","ra","ac","rf","db","sg","bh","hs","mt","ds","rg","cn","uut","uuq","uup","uuh","uus","uuo" ],
    [    0,   0,   0,"ce","pr","nd","pm","sm","eu","gd","tb","dy", "ho", "er", "tm", "yb", "lu",    0 ],
    [    0,   0,   0,"th","pa", "u","np","pu","am","cm","bk","cf", "es", "fm", "md", "no", "lr",    0 ]
];

var FULLNAME = {"h":"1 Hydrogen","he":"2 Helium","li":"3 Lithium","be":"4 Beryllium","b":"5 Boron","c":"6 Carbon","n":"7 Nitrogen","o":"8 Oxygen","f":"9 Fluorine","ne":"10 Neon","na":"11 Sodium","mg":"12 Magnesium","al":"13 Aluminium","si":"14 Silicon","p":"15 Phosphorus","s":"16 Sulfur","cl":"17 Chlorine","ar":"18 Argon","k":"19 Potassium","ca":"20 Calcium","sc":"21 Scandium","ti":"22 Titanium","v":"23 Vanadium","cr":"24 Chromium","mn":"25 Manganese","fe":"26 Iron","co":"27 Cobalt","ni":"28 Nickel","cu":"29 Copper","zn":"30 Zinc","ga":"31 Gallium","ge":"32 Germanium","as":"33 Arsenic","se":"34 Selenium","br":"35 Bromine","kr":"36 Krypton","rb":"37 Rubidium","sr":"38 Strontium","y":"39 Yttrium","zr":"40 Zirconium","nb":"41 Niobium","mo":"42 Molybdenum","tc":"43 Technetium","ru":"44 Ruthenium","rh":"45 Rhodium","pd":"46 Palladium","ag":"47 Silver","cd":"48 Cadmium","in":"49 Indium","sn":"50 Tin","sb":"51 Antimony","te":"52 Tellurium","i":"53 Iodine","xe":"54 Xenon","cs":"55 Caesium","ba":"56 Barium","la":"57 Lanthanum","ce":"58 Cerium","pr":"59 Praseodymium","nd":"60 Neodymium","pm":"61 Promethium","sm":"62 Samarium","eu":"63 Europium","gd":"64 Gadolinium","tb":"65 Terbium","dy":"66 Dysprosium","ho":"67 Holmium","er":"68 Erbium","tm":"69 Thulium","yb":"70 Ytterbium","lu":"71 Lutetium","hf":"72 Hafnium","ta":"73 Tantalum","w":"74 Tungsten","re":"75 Rhenium","os":"76 Osmium","ir":"77 Iridium","pt":"78 Platinum","au":"79 Gold","hg":"80 Mercury","tl":"81 Thallium","pb":"82 Lead","bi":"83 Bismuth","po":"84 Polonium","at":"85 Astatine","rn":"86 Radon","fr":"87 Francium","ra":"88 Radium","ac":"89 Actinium","th":"90 Thorium","pa":"91 Protactinium","u":"92 Uranium","np":"93 Neptunium","pu":"94 Plutonium","am":"95 Americium","cm":"96 Curium","bk":"97 Berkelium","cf":"98 Californium","es":"99 Einsteinium","fm":"100 Fermium","md":"101 Mendelevium","no":"102 Nobelium","lr":"103 Lawrencium","rf":"104 Rutherfordium","db":"105 Dubnium","sg":"106 Seaborgium","bh":"107 Bohrium","hs":"108 Hassium","mt":"109 Meitnerium","ds":"110 Darmstadtium","rg":"111 Roentgenium","cn":"112 Copernicium","uut":"113 Ununtrium","uuq":"114 Ununquadium","uup":"115 Ununpentium","uuh":"116 Ununhexium","uus":"117 Ununseptium","uuo":"118 Ununoctium"};

var ZFRONT = 200; // Z-plane front
var ZBACK = -200; // Z-plane back
var ZRANGE = ZFRONT - ZBACK;

var TAP_DELTA = 0.075;
var TAP_ON = 0;

// global datastore
var M_COUNTER = 0;

// animation frame wait time
var TIMEOUT = 1;

// stats
var MOVING_AVG_N = 10;     // how many samples for moving average FPS
var STATS_REFRESH_N = 10;   // how often show FPS stats

// constants
var LOG_NORMAL = 0, LOG_ERROR = 1;

// browsers
var USER_AGENT = navigator.userAgent.toLowerCase();

var IS_OPERA   = USER_AGENT.indexOf("opera") != -1;
var IS_FIREFOX = USER_AGENT.indexOf("firefox") != -1 || USER_AGENT.indexOf("minefield") != -1;
var IS_CHROME  = USER_AGENT.indexOf("chrom") != -1;
var IS_SAFARI  = USER_AGENT.indexOf("safari") != -1 && !IS_CHROME;

var IS_MAC     = USER_AGENT.indexOf("os x") != -1 || USER_AGENT.indexOf("macintosh") != -1;

/*****************************************************************************************/
// Parser
/*****************************************************************************************/
function parse_sdf(text) {
    var atoms = [];
    var bonds = [];
    var histogram = {};

    try {
        var lines = text.split("\n");

        var natoms = parseInt(lines[3].substr(0,3));
        var nbonds = parseInt(lines[3].substr(3,3));

        var x,y,z,e, first_bond = 4+natoms;
        for(var i=4; i<first_bond; ++i) {
            x = parseFloat(lines[i].substr(0,10));
            y = parseFloat(lines[i].substr(10,10));
            z = parseFloat(lines[i].substr(20,10));
            e = trim(lines[i].substr(30,5)).toLowerCase();
            // flip z to have correct handiness
			//atoms.push([x,y,z, CPK[e], capitalize(e)]);
			atoms.push([x,y,-z, CPK[e], capitalize(e)]);
            if(histogram[e]==undefined) histogram[e] = 1;
            else histogram[e] += 1;
        }

        var start, end, n;
        for(var i=first_bond; i<(first_bond+nbonds); ++i) {
            start = parseInt(lines[i].substr(0,3));
            end   = parseInt(lines[i].substr(3,3));
            n     = parseInt(lines[i].substr(6,3));
            bonds.push([start-1, end-1, n]);
        }

        return {"ok":1, "atoms":atoms, "bonds":bonds, "histogram":histogram };
    }
    catch(e) {
        log(e, LOG_ERROR);
        return {"ok":0 };
    }
}

function parse_pdb(text) {
    function hash(s,e) {
        return "s"+Math.min(s,e)+"e"+Math.max(s,e);
    }

    function parse_bond(start,length) {
        var eatom = parseInt(lines[i].substr(start,length));
        if(eatom) {
            var h = hash(satom,eatom);
            if(bhash[h]==undefined) {
                bonds.push([satom-1, eatom-1, 1]);
                bhash[h] = bonds.length-1;
            }
            else {
                // doesn't really work as almost all PDBs
                // have just normal bonds appearing multiple
                // times instead of being double/triple bonds
                //bonds[bhash[h]][2] += 1;
            }
        }
    }

    var atoms = [];
    var bonds = [];
    var histogram = {};

    var bhash = {};

    var lines = text.split("\n");
    var x,y,z,e;
    for(var i=0; i<lines.length; ++i) {
        if(lines[i].substr(0,4)=="ATOM" || lines[i].substr(0,6)=="HETATM" ) {
            x = parseFloat(lines[i].substr(30,7));
            y = parseFloat(lines[i].substr(38,7));
            z = parseFloat(lines[i].substr(46,7));
            e = trim(lines[i].substr(76,2)).toLowerCase();
            if(e=="") e = trim(lines[i].substr(12,2)).toLowerCase();
            // flip z to have correct handiness
			//atoms.push([x,y,z, CPK[e], capitalize(e)]);
			atoms.push([x,y,-z, CPK[e], capitalize(e)]);
            if(histogram[e]==undefined) histogram[e] = 1;
            else histogram[e] += 1;
        }
        else if(lines[i].substr(0,6)=="CONECT") {
            var satom = parseInt(lines[i].substr(6,5));
            parse_bond(11,5);
            parse_bond(16,5);
            parse_bond(21,5);
            parse_bond(26,5);
        }
    }

    return {"ok":1, "atoms":atoms, "bonds":bonds, "histogram":histogram };
}

/*****************************************************************************************/
// Geometry
/*****************************************************************************************/
function clamp(val, min, max) {
    if(val<min) return min;
    if(val>max) return max;
    return val;
}

function bbox(points) {
    if(points.length>0) {
        var minx = points[0][0];
        var maxx = minx;
        var miny = points[0][1];
        var maxy = miny;
        var minz = points[0][2];
        var maxz = minz;
        for(var i=1; i<points.length; ++i) {
            if(points[i][0]<minx) minx = points[i][0];
            else if(points[i][0]>maxx) maxx = points[i][0];

            if(points[i][1]<miny) miny = points[i][1];
            else if(points[i][1]>maxy) maxy = points[i][1];

            if(points[i][2]<minz) minz = points[i][2];
            else if(points[i][2]>maxz) maxz = points[i][2];
        }
        return {'x':[minx,maxx], 'y':[miny,maxy], 'z':[minz,maxz] };
    }
    else
        return {'x':[0,0], 'y':[0,0], 'z':[0,0] };
}

function translate(points, dx,dy,dz) {
    for(var i=0; i<points.length; ++i) {
        points[i][0] += dx;
        points[i][1] += dy;
        points[i][2] += dz;
    }
}

function rescale(points, scale) {
    for(var i=0; i<points.length; ++i) {
        points[i][0] *= scale;
        points[i][1] *= scale;
        points[i][2] *= scale;
    }
}

function rotate(src, dst, angles) {
    var size = src.length;

    // X-axis
    var cos = Math.cos(angles[0]);
    var sin = Math.sin(angles[0]);

    for(var i=0; i<size; ++i) {
        dst[i][0] = src[i][0];
        dst[i][1] = src[i][1] * cos - src[i][2] * sin;
        dst[i][2] = src[i][1] * sin + src[i][2] * cos;
    }

    // Y-axis
    cos = Math.cos(angles[1]);
    sin = Math.sin(angles[1]);

    var a,b;
    for(i=0; i<size; ++i) {
        a = dst[i][0] * cos - dst[i][2] * sin;
        b = dst[i][0] * sin + dst[i][2] * cos;
        dst[i][0] = a;
        dst[i][2] = b;
    }

    // Z-axis
    cos = Math.cos(angles[2]);
    sin = Math.sin(angles[2]);

    for(i=0; i<size; ++i) {
        a = dst[i][0] * cos - dst[i][1] * sin;
        b = dst[i][0] * sin + dst[i][1] * cos;
        dst[i][0] = a;
        dst[i][1] = b;
    }
}

/*****************************************************************************************/
// Canvas rendering
/*****************************************************************************************/
function refresh_background(ctx, width, height, color) {
    ctx.fillStyle = color;
    ctx.fillRect(0, 0, width, height);
}

function point3d(ctx, x,y,z, radius) {
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, 6.28, 0);
    ctx.fill();
}

function str3d(ctx, str, x,y,z, radius) {
    ctx.font = (2.4*radius)+"px Arial";
    ctx.fillText(str, x-0.7*radius, y+0.7*radius);
}

function line3d(ctx, x1,y1,z1, x2,y2,z2, scolor,ecolor, bonds_gradient) {
    ctx.beginPath();
    ctx.moveTo(x1,y1);
    ctx.lineTo(x2,y2);
    ctx.closePath();

    if(bonds_gradient) {
        var lingrad = ctx.createLinearGradient(x1,y1,x2,y2);
        lingrad.addColorStop(0.0, scolor);
        lingrad.addColorStop(1.0, ecolor);
        ctx.strokeStyle = lingrad;
    }
    else {
        ctx.strokeStyle = scolor;
    }
    ctx.stroke();
}

function render_molecule(ctx, cx, cy, width, height, points, lines, draw_atoms, draw_bonds,
                         bonds_color, bonds_autowidth, bonds_gradient,
                         atom_colors, atom_radius, flat_shading, sphere_shading, letters_shading,
                         zback, zrange) {
    var r,g,b, r1,g1,b1, r2,g2,b2;
    var c1,c2, col1,col2, elcol, x,y,z, d, z1,z2;
    var s,dx,dy,dz,l;
    var grad;
    var half, size;
    var dxx, dyy, dd, fi;

    var atom_margin_coef = 0.85;
    if(letters_shading) atom_margin_coef = 1.2;

    // Z-sort atoms and bonds
    var elements = [];
    for(var i=0; i<points.length; ++i)
        elements.push([points[i][2], 0, i]);
    for(i=0; i<lines.length; ++i) {
        var midpoint = 0.5*(points[lines[i][0]][2] + points[lines[i][1]][2]);
        elements.push([midpoint, 1, i, lines[i][2]]);
    }
    elements.sort(function(a,b) { return a[0]-b[0]; } );

    for(i=0; i<elements.length; ++i) {

        // render atoms
        if(draw_atoms && elements[i][1]==0) {
            x = cx + points[elements[i][2]][0];
            y = cy + points[elements[i][2]][1];

            // cull points outside of canvas
            if(x>=-atom_radius && x<=width+atom_radius && y>=-atom_radius && y<=height+atom_radius) {

                // set color depending on depth
                z = points[elements[i][2]][2];
                col1 = clamp(Math.round(255*(z-zback)/zrange), 0, 255);

                if(atom_colors) {
                    elcol = points[elements[i][2]][3];
                    r1 = elcol[0];
                    g1 = elcol[1];
                    b1 = elcol[2];

                }
                else {
                    r1 = ATOMS_COLOR[0];
                    g1 = ATOMS_COLOR[1];
                    b1 = ATOMS_COLOR[2];
                }

                r2 = BACKGROUND_COLOR[0];
                g2 = BACKGROUND_COLOR[1];
                b2 = BACKGROUND_COLOR[2];

                d = (255-col1)/255;

                r = Math.round(r1+d*(r2-r1));
                g = Math.round(g1+d*(g2-g1));
                b = Math.round(b1+d*(b2-b1));

                if(letters_shading) {
                    ctx.fillStyle ='rgb('+r+','+g+','+b+')';
                    str3d(ctx, points[elements[i][2]][4], x, y, z, (0.5+col1/512)*atom_radius);
                }
                else if(flat_shading) {
                    ctx.fillStyle ='rgb('+r+','+g+','+b+')';
                    point3d(ctx, x, y, z, (0.5+col1/512)*atom_radius);
                }
                else if(sphere_shading) {
                    half = (0.5+col1/512)*atom_radius;
                    size = 2*half;
                    grad = ctx.createRadialGradient(x-half+size*0.34,y-half+size*0.7,size*0.15, x,y,half);
                    grad.addColorStop(0, 'rgb('+r+','+g+','+b+')');
                    grad.addColorStop(0.95, 'rgba(0,0,0,1.0)');
                    grad.addColorStop(1, 'rgba(0,0,0,0.0)');

                    ctx.fillStyle = grad;

                    // Safari doesn't handle well transparent part of the radial gradient
                    // so it has to be clipped by circle shape.
                    // Other browsers can just draw rectangle which is faster.
                    if(IS_SAFARI)
                        point3d(ctx, x, y, z, (0.5+col1/512)*atom_radius);
                    else
                        ctx.fillRect(x-half, y-half, size,size);
                }
                else {
                    dx = atom_radius*1.5;

                    // must invert circle and colors order for Opera
                    if(IS_OPERA) {
                        grad = ctx.createRadialGradient(x+dx, y+dx, 1, x,y,atom_radius);
                        grad.addColorStop(0.0, 'rgb(0,0,0)');
                        grad.addColorStop(1.0, 'rgb('+r+','+g+','+b+')');
                    }

                    else {
                        grad = ctx.createRadialGradient(x,y,atom_radius, x+dx, y+dx, 1);
                        grad.addColorStop(0.0, 'rgb('+r+','+g+','+b+')');
                        grad.addColorStop(1.0, 'rgb(0,0,0)');
                    }

                    ctx.fillStyle = grad;
                    point3d(ctx, x, y, z, (0.5+col1/512)*atom_radius);
                }
            }
        }

        // render bonds
        else if(draw_bonds && elements[i][1]==1) {
            x1 = points[lines[elements[i][2]][0]][0] + cx;
            y1 = points[lines[elements[i][2]][0]][1] + cy;
            z1 = points[lines[elements[i][2]][0]][2];
            x2 = points[lines[elements[i][2]][1]][0] + cx;
            y2 = points[lines[elements[i][2]][1]][1] + cy;
            z2 = points[lines[elements[i][2]][1]][2];
            col1 = clamp(Math.round(255*(z1-zback)/zrange), 0, 255);
            col2 = clamp(Math.round(255*(z2-zback)/zrange), 0, 255);

            if(bonds_autowidth)
                ctx.lineWidth = 0.1+col1/120.0;

            dx = x2 - x1;
            dy = y2 - y1;
            dz = z2 - z1;
            l = Math.sqrt(dx*dx+dy*dy+dz*dz);

            if(draw_atoms) {
                s = atom_margin_coef*atom_radius*dx/l;
                x1 += s;
                x2 -= s;

                s = atom_margin_coef*atom_radius*dy/l;
                y1 += s;
                y2 -= s;
            }

            r = bonds_color[0];
            g = bonds_color[1];
            b = bonds_color[2];
            c1 = "rgba("+r+","+g+","+b+","+col1/255+")";
            c2 = "rgba("+r+","+g+","+b+","+col2/255+")";

            // first bond
            line3d(ctx, x1,y1,z1, x2,y2,z2,  c1,c2, bonds_gradient);

            // second and third bonds
            if(elements[i][3]>1) {
                dd = 0.35*atom_radius;
                fi = Math.atan2(dy, dx)-1.5;
                dxx = dd*Math.cos(fi);
                dyy = dd*Math.sin(fi);

                if(!draw_atoms) {
                    s = 0.6*atom_radius*dx/l;
                    x1 += s;
                    x2 -= s;

                    s = 0.6*atom_radius*dy/l;
                    y1 += s;
                    y2 -= s;
                }

                if(elements[i][3]==2) {
                    line3d(ctx, x1+dxx,y1+dyy,z1, x2+dxx,y2+dyy,z2, c1,c2, bonds_gradient);
                }
                else if(elements[i][3]==2) {
                    line3d(ctx, x1+dxx,y1+dyy,z1, x2+dxx,y2+dyy,z2, c1,c2, bonds_gradient);
                    line3d(ctx, x1-dxx,y1-dyy,z1, x2-dxx,y2-dyy,z2, c1,c2, bonds_gradient);
                }
            }
        }
    }

}

/*****************************************************************************************/
// Molecule
/*****************************************************************************************/
function Molecule(pars) {
    this.mdiv = pars.mdiv;

    this.canvas = document.getElementById(pars.canvas_id);
    this.ctx = this.canvas.getContext("2d");
    this.width = this.canvas.width;
    this.height = this.canvas.height;
    this.cx = this.width/2;
    this.cy = this.height/2;

    this.angles = [0, 0, 0];
    this.points = [];

    this.atoms = [];
    this.bonds = [];

    this.draw_atoms = pars.draw_atoms;
    this.draw_bonds = pars.draw_bonds;

    this.bonds_autowidth = pars.bonds_autowidth;
    this.bonds_gradient = pars.bonds_gradient;

    this.flat_shading = pars.flat_shading;
    this.sphere_shading = pars.sphere_shading;
    this.letters_shading = 0;

    this.atom_colors = pars.atom_colors;

    this.bonds_color = BONDS_COLOR;

    this.autorotation = [0,0,0]; // Autorotation around x,y,z axes
    this.timeout_id = -1;

    this.drag_on = 0;
    this.drag_startx = 0;
    this.drag_starty = 0;

    this.show_stats = 0;
    this.stats_cache = "";
    this.stats_refresh = 0;

    this.stats_count = 0;
    this.stats_total = 0;

    this.stats_buf = new Array(MOVING_AVG_N);
    this.stats_index = 0;

    this.time_start = 0;

    this.atom_radius = pars.atom_radius;
    this.zback = ZBACK;
    this.zrange = ZRANGE;

    this.destruct = function() {
        if(this.timeout_id>=0)
            clearTimeout(this.timeout_id);
    }

    this.render = function() {
        refresh_background(this.ctx, this.width, this.height, BACKGROUND_COLOR_RGB);
        render_molecule(this.ctx, this.cx, this.cy, this.width, this.height,
                        this.points, this.bonds, this.draw_atoms, this.draw_bonds,
                        this.bonds_color, this.bonds_autowidth, this.bonds_gradient,
                        this.atom_colors, this.atom_radius,
                        this.flat_shading, this.sphere_shading, this.letters_shading,
                        this.zback, this.zrange);
    }

    this.reset_model = function() {
        var bb = bbox(this.atoms);

        // center model (middle of bounding box)
        var cx = bb.x[0]+(bb.x[1]-bb.x[0])/2.0;
        var cy = bb.y[0]+(bb.y[1]-bb.y[0])/2.0;
        var cz = bb.z[0]+(bb.z[1]-bb.z[0])/2.0;
        translate(this.atoms, -cx,-cy,-cz);

        // scale model (fit in canvas)
        var max_model = Math.max(bb.x[1]-bb.x[0], bb.y[1]-bb.y[0]);
        var min_canvas = Math.min(this.width, this.height);
        var scale = 0.8*min_canvas/max_model;
        rescale(this.atoms, scale);

        // copy original atom locations
        this.points = array2d(this.atoms.length, 4);
        for(var i=0; i<this.atoms.length; i++) {
            this.points[i][0] = this.atoms[i][0];
            this.points[i][1] = this.atoms[i][1];
            this.points[i][2] = this.atoms[i][2];
            this.points[i][3] = this.atoms[i][3];
            this.points[i][4] = this.atoms[i][4];
        }
        rotate(this.atoms, this.points, this.angles);

        this.stats_cache = generate_stats(this.atoms.length, this.bonds.length);
        this.mdiv.children(".dragheader").children(".stats").text( this.stats_cache );

        //this.highlight_elements();
    }

    this.highlight_elements = function() {
        var elements = [];
        for(var i in this.histogram) elements.push(i);
        highlight_elements(elements);
    }

    this.process_sdf = function(data) {
        var result = parse_sdf(data);
        if(result.ok) {
            this.atoms = result.atoms;
            this.bonds = result.bonds;
            this.histogram = result.histogram;

            this.reset_model();
            this.render();
        }
    }

    this.process_pdb = function(data) {
        var result = parse_pdb(data);
        if(result.ok) {
            this.atoms = result.atoms;
            this.bonds = result.bonds;
            this.histogram = result.histogram;

            this.reset_model();
            this.render();
        }
    }

    this.animate = function() {
        if(this.show_stats)
            var start = (new Date).getTime();

        if(this.autorotation[0]) this.angles[0] += 0.01;
        if(this.autorotation[1]) this.angles[1] += 0.01;
        if(this.autorotation[2]) this.angles[2] += 0.01;
        rotate(this.atoms, this.points, this.angles);
        this.render();

        if(this.show_stats) {
            // FPS based on real elapsed time
            if(this.stats_refresh == STATS_REFRESH_N) {
                var current = (new Date).getTime();
                var diff = current - this.time_start;
                var fps = 1000*this.stats_refresh/diff;

                var stats = this.stats_cache + " " + fps.toFixed(1) + " fps";
                this.mdiv.children(".dragheader").children(".stats").text( stats );

                this.stats_refresh = 0;
                this.time_start = current;
            }
            this.stats_refresh += 1;

            /*
            // FPS computed from pure rendering & computation time
            var diff = (new Date).getTime() - start;

            // moving average FPS stats
            // keep track of total sum and count
            this.stats_total += diff;
            this.stats_count += 1;
            this.stats_refresh += 1;

            // store last N values
            this.stats_buf[this.stats_index] = diff;
            this.stats_index += 1;
            if(this.stats_index == MOVING_AVG_N)
                this.stats_index = 0;

            // render stats only every X frames
            if(this.stats_count == MOVING_AVG_N && this.stats_refresh == STATS_REFRESH_N) {
                var stats_avg = this.stats_total/this.stats_count;
                var stats = this.stats_cache + " " + (1000.0/stats_avg).toFixed(1) + " fps";
                this.mdiv.children(".dragheader").children(".stats").text( stats );
            }
            if(this.stats_refresh == STATS_REFRESH_N)
                this.stats_refresh = 0;

            // subtract oldest value
            if(this.stats_count == MOVING_AVG_N) {
                this.stats_total -= this.stats_buf[this.stats_index];
                this.stats_count -= 1;
            }
            */

            /*
            // simple FPS stats (blocks of N values)
            if(this.stats_count<10) {
                this.stats_total += diff;
                this.stats_count++;
            }
            else {
                var stats_avg = this.stats_total/this.stats_count;
                var stats = this.stats_cache + " " + (1000.0/stats_avg).toFixed(1) + " fps";
                this.mdiv.children(".dragheader").children(".stats").text( stats );
                this.stats_count = 0;
                this.stats_total = 0;
            }
            */
        }
    }

    this.toggle_autorotation = function(axis, what) {
        this.autorotation[axis] =! this.autorotation[axis];

        if(this.timeout_id>=0)
            clearTimeout(this.timeout_id);

        var auto_on = this.autorotation[0] || this.autorotation[1] || this.autorotation[2];
        if(auto_on)
            this.timeout_id = setInterval(function() { what.animate() }, TIMEOUT);
    }

    this.wheel = function(e,d) {
        var scale = 1+0.05*d;
        this.atom_radius *= scale;
        this.zback *= scale;
        this.zrange *= scale;
        rescale(this.atoms, scale);
        rotate(this.atoms, this.points, this.angles);
        this.render();
        e.stopPropagation();
        e.preventDefault();
    }

    this.drag_start = function(e) {
        this.drag_on = 1;
        this.drag_startx = e.pageX;
        this.drag_starty = e.pageY;
        e.stopPropagation();
    }

    this.drag_end = function(e) {
        this.drag_on = 0;

        this.angles = [0, 0, 0];
        for(var i=0; i<this.atoms.length;i++) {
            var c = this.points[i];
            this.atoms[i] = [c[0], c[1], c[2]];
        }
        e.stopPropagation();
    }

    this.drag = function(e) {
        if(this.drag_on) {
            var deltaX = this.drag_startx - e.pageX;
            var deltaY = this.drag_starty - e.pageY;

            this.drag_startx = e.pageX;
            this.drag_starty = e.pageY;

            this.angles[1] += 0.01*deltaX;
            this.angles[0] += 0.01*deltaY;

            rotate(this.atoms, this.points, this.angles);
            this.render();
        }
        e.stopPropagation();
    }

    this.click_left = function() {
        this.angles[1] += TAP_DELTA;
        rotate(this.atoms, this.points, this.angles);
        this.render();
    }

    this.click_right = function() {
        this.angles[1] += -TAP_DELTA;
        rotate(this.atoms, this.points, this.angles);
        this.render();
    }

    this.click_up = function() {
        this.angles[0] += TAP_DELTA;
        rotate(this.atoms, this.points, this.angles);
        this.render();
    }

    this.click_down = function() {
        this.angles[0] += -TAP_DELTA;
        rotate(this.atoms, this.points, this.angles);
        this.render();
    }

    this.toggle_atoms = function() {
        this.draw_atoms = !this.draw_atoms;
        this.render();
    }

    this.toggle_bonds = function() {
        this.draw_bonds = !this.draw_bonds;
        this.render();
    }

    this.toggle_autorotation_x = function() {
        this.toggle_autorotation(0, this);
    }

    this.toggle_autorotation_y = function() {
        this.toggle_autorotation(1, this);
    }

    this.toggle_autorotation_z = function() {
        this.toggle_autorotation(2, this);
    }

    this.toggle_width = function() {
        this.bonds_autowidth = !this.bonds_autowidth;
        this.render();
    }

    this.toggle_gradient = function() {
        this.bonds_gradient = !this.bonds_gradient;
        this.render();
    }

    this.toggle_color = function() {
        this.atom_colors = !this.atom_colors;
        this.render();
    }

    this.toggle_flat = function() {
        this.flat_shading = !this.flat_shading;
        this.render();
    }

    this.toggle_sphere = function() {
        this.sphere_shading = !this.sphere_shading;
        this.render();
    }

    this.toggle_letters = function() {
        this.letters_shading = !this.letters_shading;
        this.render();
    }
}

/*****************************************************************************************/
// Ajax
/*****************************************************************************************/
function ajax_load(url, type, callback) {
    $.ajax({
        url: url,
        dataType: type,
        success: function(data) {
            if(data) {
                if(callback) callback(data);
            }
        },
        error: function(XMLHttpRequest, textStatus, errorThrown) {
            var message = "Error ["+textStatus+"]["+errorThrown+"]";
            log(message, LOG_ERROR);
        }
    });
}

/*****************************************************************************************/
// Helpers
/*****************************************************************************************/
function log(message, type) {
    $("#log").css("display","block");
    if(type && type == LOG_ERROR)
        message = "<span class='error'>"+message+"</span>";
    $("#log").append(message+"<br/>");
}

function trim(text) {
    return text.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
}

function capitalize(text) {
    return text.charAt(0).toUpperCase() + text.substr(1).toLowerCase();
}

function escape(txt) {
    return txt.replace(/ /g,"_");
}

function unescape(txt) {
    return txt.replace(/_/g," ");
}

function pretty_number(x) {
    var delimiter = "'";
    var strx = x.toString();
    var pretty = "";
    for(var i=strx.length-1; i>=0; i-- ) {
        if((strx.length-1-i) % 3 == 0 && (strx.length-1-i) > 0)
            pretty = delimiter+pretty;
        pretty = strx[i]+pretty;
    }
    return pretty;
}

function pluralize(count) {
    if(count==1) return "";
    return "s";
}

function array2d(rows, cols) {
    var r = new Array(rows);
    for(var c=0; c<r.length; ++c)
        r[c] = new Array(cols);
    return r;
}

function bind(scope, fn) {
    return function () {
        fn.apply(scope, arguments);
    };
}

function guess_type(text) {
    var lines = text.split("\n");
    for(var i=0; i<lines.length; ++i) {
        if(lines[i].substr(0,4)=="ATOM" || lines[i].substr(0,6)=="HETATM" ) {
            return "pdb"
        }
    }
    return "sdf";
}

/*****************************************************************************************/
// HTML rendering
/*****************************************************************************************/
function generate_stats(natoms, nbonds) {
    return pretty_number(natoms)+" atom"+pluralize(natoms)+" "+pretty_number(nbonds)+" bond"+pluralize(nbonds);
}

function create_permalink(txt) {
    return '<a href="#'+escape(txt)+'">'+txt+'</a>';
}

function molecule_div(i) {
    var tmp = '';
    tmp += '<div class="molecule shadow">';
    tmp += '<div class="dragheader">Drag me</div>';
    tmp += '<canvas width="400" height="400" id="screen'+i+'"></canvas>';
    tmp += '<br class="clr"/>';

    tmp += '<div class="b_left touch"> </div>';
    tmp += '<div class="b_right touch"> </div>';
    tmp += '<div class="b_up touch"> </div>';
    tmp += '<div class="b_down touch"> </div>';

    tmp += '<div class="button b_x">X</div>';
    tmp += '<div class="button b_y">Y</div>';
    tmp += '<div class="button b_z">Z</div>';
    tmp += '<div class="button b_atoms" title="Toggle atoms">Atm</div>';
    tmp += '<div class="button b_bonds" title="Toggle bonds">Bnd</div>';
    tmp += '<div class="button b_width" title="Bonds width">BW</div>';
    tmp += '<div class="button b_grad" title="Bonds gradient">BG</div>';
    tmp += '<div class="button b_col" title="Atoms colors">ACol</div>';
    tmp += '<div class="button b_flt" title="Atoms flat shading">AFlt</div>';
    tmp += '<div class="button b_sph" title="Atoms sphere shading">ASph</div>';
    tmp += '<div class="button b_ltr" title="Atoms letters shading">ALtr</div>';
    tmp += '<div class="button b_png"><a href="#">PNG</a></div>';
    tmp += '<div class="button b_stats">Stats</div>';
    tmp += '</div>';
    return tmp;
}

function add_molecule_canvas(label, name, text, target) {
    var i = M_COUNTER;
    M_COUNTER++;

    var mdiv = $(molecule_div(i));
    //$("#molecules").append(mdiv);
    target.append(mdiv);
    if(TAP_ON) mdiv.children(".touch").css("display","block");

    var pars = PARS[get_type(name)];
    pars['canvas_id'] = "screen"+i;
    pars['mdiv'] = mdiv;

    // sorry for some reason Safari Win performance sucks incredibly for radial gradients
    // also Firefox on OSX was reported to hang up
    if((IS_SAFARI && !IS_MAC) || (IS_FIREFOX && IS_MAC)) pars.flat_shading = 1;

    var canvas = $("#"+pars.canvas_id);
    var parent = canvas.parent();

    var header = '<span class="close">X</span>'+label+'<span class="stats">X atoms Y bonds Z FPS</span>';
    parent.children(".dragheader").html(header);

    var molecule = new Molecule(pars);

    var parser_map = {
    "pdb":molecule.process_pdb,
    "sdf":molecule.process_sdf,
    "mol":molecule.process_sdf
    }
    if(text) {
        var ext = guess_type(text);
        bind(molecule, parser_map[ext])(text);
    }
    else {
        var url = "data/"+name;
        var ext = name.substr(-3);
        ajax_load(url, "txt", bind(molecule, parser_map[ext]) );
    }


    canvas.mousedown( bind(molecule, molecule.drag_start) );
    $(document).mouseup( bind(molecule, molecule.drag_end) );
    $(document).mousemove( bind(molecule, molecule.drag) );

    canvas.mousewheel( bind(molecule, molecule.wheel) );

    if(pars.draw_atoms) parent.children(".b_atoms").addClass("on");
    if(pars.draw_bonds) parent.children(".b_bonds").addClass("on");

    if(pars.bonds_autowidth) parent.children(".b_width").addClass("on");
    if(pars.bonds_gradient)  parent.children(".b_grad").addClass("on");
    if(pars.atom_colors)     parent.children(".b_col").addClass("on");
    if(pars.flat_shading)    parent.children(".b_flt").addClass("on");
    if(pars.sphere_shading)  parent.children(".b_sph").addClass("on");

    if(pars.auto_x) {
        parent.children(".b_x").addClass("on");
        molecule.toggle_autorotation(0, molecule);
    }
    if(pars.auto_y) {
        parent.children(".b_y").addClass("on");
        molecule.toggle_autorotation(1, molecule);
    }
    if(pars.auto_z) {
        parent.children(".b_z").addClass("on");
        molecule.toggle_autorotation(2, molecule);
    }

    //$(".molecule").draggable();

    parent.children(".b_bonds").click( bind(molecule, molecule.toggle_bonds ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_atoms").click( bind(molecule, molecule.toggle_atoms ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });

    parent.children(".b_x").click( bind(molecule, molecule.toggle_autorotation_x ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_y").click( bind(molecule, molecule.toggle_autorotation_y ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_z").click( bind(molecule, molecule.toggle_autorotation_z ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });

    parent.children(".b_width").click( bind(molecule, molecule.toggle_width ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_grad").click( bind(molecule, molecule.toggle_gradient ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });

    parent.children(".b_col").click( bind(molecule, molecule.toggle_color ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_flt").click( bind(molecule, molecule.toggle_flat ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_sph").click( bind(molecule, molecule.toggle_sphere ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });
    parent.children(".b_ltr").click( bind(molecule, molecule.toggle_letters ) ).toggle( function() { $(this).toggleClass("on") }, function() { $(this).toggleClass("on") });

    parent.children(".b_stats").toggle( function() { parent.children(".dragheader").children(".stats").css("display", "inline"); $(this).toggleClass("on"); molecule.show_stats = 1; },
                                        function() { parent.children(".dragheader").children(".stats").css("display", "none");   $(this).toggleClass("on"); molecule.show_stats = 0; } );

    //parent.children(".b_png").click( function() { window.location.href = canvas.get(0).toDataURL("image/png").replace("image/png","image/octet-stream") } );
    parent.children(".b_png").children("a").click( function() { this.href = canvas.get(0).toDataURL("image/png") } );

    parent.children(".b_left").click( bind(molecule, molecule.click_left) );
    parent.children(".b_right").click( bind(molecule, molecule.click_right) );
    parent.children(".b_up").click( bind(molecule, molecule.click_up) );
    parent.children(".b_down").click( bind(molecule, molecule.click_down) );

    parent.children(".dragheader").children(".close").click( function() {
        parent.remove();

        molecule.destruct();
        delete molecule;
    });

    /* Doesn't work, somewhere data gets corrupted
    $(".molecule").resizable({
        resize: function(event, ui) {
            molecule.width = $(this).width();
            molecule.height = $(this).height()-20;
            canvas.attr('width', molecule.width).attr('height', molecule.height);

            molecule.cx = molecule.width/2;
            molecule.cy = molecule.height/2;

            molecule.reset_model();
            molecule.render();
        },
        start: function(event, ui) {
            //molecule.pause_rotation(molecule);
        },
        stop: function(event, ui) {
            //molecule.restore_rotation(molecule);
        }
    });
    */
}

/*****************************************************************************************/
// Periodic table
/*****************************************************************************************/
function create_periodic_table(parent_id) {
    var parent =$("#"+parent_id);
    for(var r=0; r<MENDELEEV.length; ++r) {
        for(var c=0; c<MENDELEEV[r].length; ++c) {
            var element = MENDELEEV[r][c];
            if(element) {
                var color = CPK[element];
                var name = capitalize(element);
                var div = $("<div class='mm' id='m_"+element+"' title='"+FULLNAME[element]+"' style='color:rgb("+color[0]+","+color[1]+","+color[2]+")'>"+name+"</div>");
                parent.append(div);
                div.css("textShadow", "0px 0px 4px rgb("+color[0]+","+color[1]+","+color[2]+")");
            }
            else
                parent.append($("<div>&nbsp;</div>"));
        }
        parent.append($("<br class='clr'/>"));
    }
}

function highlight_elements(elements) {
    $(".mm").addClass("neutral");
    for(var i=0; i<elements.length; ++i) {
         $("#m_"+elements[i]).removeClass("neutral");
    }
}

function reset_highlights() {
    $(".mm").removeClass("neutral");
}

/*****************************************************************************************/
// Molecule menu
/*****************************************************************************************/
function create_molecule_menu(parent_id, add_func) {
    var parent =$("#"+parent_id);
    var first = 1;
    for(var i in MENU_MAP) {
        if(MENU_MAP[i])
            parent.append($("<li class='mol'>"+i+"</li>"));
        else
            if(first) {
                parent.append($("<li class='hedtop'>"+i+"</li>"));
                first = 0;
            }
            else
                parent.append($("<li class='hed'>"+i+"</li>"));
    }

    $("#"+parent_id+" li.mol").click(function() {
        var txt = $(this).text();
        add_func(create_permalink(txt), MENU_MAP[txt]);
        window.location.hash = escape(txt);
    });
}

/*****************************************************************************************/
// Parameters
/*****************************************************************************************/
function get_parameter(url) {
    var chunks = url.split("#");
    if(chunks.length==2)
        return unescape(chunks[1].toLowerCase());
    return "";
}
