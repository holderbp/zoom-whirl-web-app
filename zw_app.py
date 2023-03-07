import dash
import dash.html as html
import dash.dcc as dcc
import dash_bootstrap_components as dbc
import dash.exceptions as dex
import plotly.express as px
import plotly.graph_objects as pgo
import numpy as np
import zw_app_html as zwah
import zw_app_orbitcalc as zwoc

from dash.dependencies import Input, Output, State

#
#  This includes an example of a client-side callback for fast animation
#  that I stole from @emher's 2020-09-01 stackexchange answer here:
#
#      https://stackoverflow.com/questions/63589249
#
#  We hope to adapt it for zoom-whirl purposes.
#


#
#--- Create the app and its html
#
stylesheets = [
    {
        "href": "https://fonts.googleapis.com/css2?"
        "family=Lato:wght@400;700&display=swap",
        "rel": "stylesheet",
    }
]
app = dash.Dash(
    __name__,
    update_title=None,  # remove "Updating..." from title
    external_stylesheets =
    [
        dbc.themes.BOOTSTRAP,
        stylesheets,
    ],
    suppress_callback_exceptions=True,  # suppress initial load errors
)
#
#--- Create a server (used by heroku)
#
server = app.server

#
# The "backbone" of the site was taken from my COVID portal.  Web
# content will be placed into "page-content" below.
#
site_backbone = html.Div(
    [
        dcc.Location(id='url', refresh=False),
        html.Div(id='page-content', className="page border")
    ]
)
app.layout = site_backbone

###############################
#      parameter values       #
###############################
#
#--- physical constants
#
#   FIXME: should move to units of r_s = 2GM (graphs in r/r_s)
#
G = 1
M = 1
#
#--- orbit graph
#
default_orbit_resolution = 1000
print_out_data = False # This should usually be false... will print out each orbit
#
#--- effective potential graph
#
# Npoints, and min and max values for making graph
Npoints_r = 20000
rmin = 0.1
rmax = 3000
# max allowed ang momentum value
ell_max = float(zwah.ell_max_str)
# for choosing an energy just below Vmax
energy_frac_of_Vmax = 0.9
# starting values
default_angmom_str = "3.7" # "4.0"
default_energy_str = "-0.03445" # "-0.03"
default_ecc_str = "0.633223" # 
default_periap_str = "4.52351" #

###############################
#       helper functions      #
###############################

def makefig_orbit(orb_r, orb_phi, ra):
    #
    #--- Create the plotly figure
    #
    fig = pgo.Figure()
    # trace [0]: static dashed trace
    fig.add_trace(
        pgo.Scatterpolar(
            r = orb_r,
            theta = orb_phi,
            thetaunit = "radians",
            mode = 'lines',
            line = dict(dash = 'dash',
                        color = 'grey',
                        width = 1),
        )
    )
    # # trace [1]: moving trace (empty for now)
    fig.add_trace(
         pgo.Scatterpolar(
             r = [orb_r[0]],
             theta = [orb_phi[0]],
             thetaunit = "radians",             
             mode = 'markers',
             marker = dict(
                 color='lightslategray',
                 size=8,
                 symbol='circle'),
         )
    )
    # update the axes labels and title
    fig.update_layout(
        polar = dict(
            radialaxis = dict(
                range=[0, 1.02*ra],
                #showticklabels=False, ticks=''
            ),
            angularaxis = dict(
                #showticklabels=False, ticks=''
            )
        ),
        margin=dict(l=35, r=10, t=10, b=10),
        showlegend=False,
    )
    return fig

def create_orbit_figure(ell, E):
    # Set the values of ang mom and energy in orbitcalc module
    zwoc.ell = ell
    zwoc.E = E
    # evolve the orbit
    t, r_t, phi_t, rp, ra, ecc = zwoc.get_orbit()
    # make a figure
    fig = makefig_orbit(r_t, phi_t, ra)
    # return the figure and the data
    return fig, t, r_t, phi_t

def get_ell_E_and_check_pars(ecc, periap, version=None):
    """
    Convert (e, rp) -> (ell, E), and verify that the result is valid.
    If the result is not valid, try to change either eccentricity or
    periapsis slightly (whichever was not set by user) to make it so.
    """
    if version == 'ecc':
        # check validity:
        #     ell**2 > 0
        rp = periap
        while True:
            if ( (1 + ecc) * rp <= (3 + ecc**2) ):
                rp *= 1.05
            else:
                break
        # check validity:
        #    ell < ell_max   and   ell > sqrt(4)   and   E < 0
        while True:
            ell = (1 + ecc) * rp / np.sqrt( (1 + ecc) * rp - (3 + ecc**2) )
            E = ( (1 - ecc) / (2 * rp) * (4 - (1 + ecc) * rp)
                  / ( (1 + ecc) * rp - (3 + ecc**2) ) )
            if ( ell > ell_max ):
                # failure
                rp = None
                break
            elif ( (ell**2 < 12) | (E > 0) ):
                rp *=1.05
            else:
                break
        periap = rp
        return [ell, E, ecc, rp]
    elif version == 'periap':
        # check validity:
        #     ell**2 > 0
        e = ecc
        while True:
            if ( (1 + e) * periap <= (3 + e**2) ):
                e += (1-e)/20
            else:
                break
        # check validity:
        #    ell < ell_max   and   ell > sqrt(4)   and   E < 0
        while True:
            ell = (1 + e) * periap / np.sqrt( (1 + e) * periap - (3 + e**2) )
            E = ( (1 - e) / (2 * periap) * (4 - (1 + e) * periap)
                  / ( (1 + e) * periap - (3 + e**2) ) )
            if ( ell > ell_max ):
                # failure
                e = None
                break
            elif ( (ell**2 < 12) | (E > 0) ):
                e += (1-e)/20
            else:
                break
        return [ell, E, e, periap]

def get_ecc_periap_and_check_E(ell, E):
    zwoc.G = G
    zwoc.M = M
    zwoc.ell = ell
    zwoc.E = E
    # find the positions two circular orbits (inner-unstable, outer-stable)
    r_inner, r_outer = zwoc.get_Veff_maxmin_r_values()
    Vmin = zwoc.Veff(r_outer)
    Vmax = zwoc.Veff(r_inner)
    deltaV = Vmax-Vmin
    #
    #--- check the energy value to make sure it is valid
    #
    if E < Vmin:
        # just set to Vmin (stable circular orbit)
        E = Vmin
        zwoc.E = E
    elif E > Vmax:
        # forbid plunging orbits and hyperbolic orbits...
        # set to just below peak or just below zero
        #
        # FIXME: allow hyperbolic/plunging?
        #
        E1 = energy_frac_of_Vmax*deltaV + Vmin
        E2 = energy_frac_of_Vmax*(0 - Vmin) + Vmin
        E = min(E1, E2)
        zwoc.E = E
    # get periapsis and apoapsis and calculate eccentricity
    rp, ra = zwoc.get_peri_and_apoapsis()
    ecc = zwoc.get_eccentricity(rp, ra)
    return [ecc, rp, E]

def get_number_from_string(valstr):
    try:
        val = float(valstr)
        return val
    except ValueError:
        return None

def get_all_values_from_strings(angmom_str, energy_str, ecc_str, periap_str, checkvalid=False):
    ell = get_number_from_string(angmom_str)
    E = get_number_from_string(energy_str)
    ecc = get_number_from_string(ecc_str)
    periap = get_number_from_string(periap_str)
    # check validity
    if checkvalid:
        if ( (ell**2 < 12*(G*M)**2) | (ell > ell_max) ):
            ell = None
        if (E > 0):
            E = None
        if ( (ecc < 0) | (ecc >= 1) ):
            ecc = None
        if (periap < 4*G*M):
            periap = None
    return [ell, E, ecc, periap]

def makefig_effpot(r, V, VN, E, win_rmin, win_rmax, win_Vmin, win_Vmax):
    # make the plotly figure object
    fig = pgo.Figure()
    # effective potential graph
    fig.add_trace(
        pgo.Scatter(x=r, y=V, mode='lines')
    )
    # constant E line
    fig.add_trace(
        pgo.Scatter(x=r, y=E*np.ones(len(r)), mode='lines')
    )
    #
    # Newtonian potential  (this doesn't work
    #fig.add_trace(
    #pgo.Scatter(x=r, y=VN)#, mode='lines')
    #)
    #
    # moving point on E line
    #  (empty.. will be filled by orbit eval of E)
    fig.add_trace(
        pgo.Scatter(x=[], y = [],
                    marker = dict(
                        color='lightslategray',
                        size=8,
                        symbol='circle'),
                    )
    )
    fig.update_layout(
        #title_text = 'Effective Potential',
        xaxis = dict(range = [win_rmin, win_rmax],
                     title = {'text': 'Radius/M'}),
        yaxis = dict(range = [win_Vmin, win_Vmax],
                     title = {'text': 'Veff'}),
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,        
    )
    return fig

def makefig_gw(t, H, Htype):
    # make the plotly figure object
    fig = pgo.Figure()
    # graph of H (plus or cross)
    fig.add_trace(
        pgo.Scatter(
            x=t, y=H,
            mode = 'lines',
            line = dict(dash = 'dash',
                        color = 'grey',
                        width = 1),
            )
    )
    # moving point on graph
    #  (empty.. will be filled by orbit eval of E)
    fig.add_trace(
        pgo.Scatter(x=[], y = [],
                    marker = dict(
                        color='lightslategray',
                        size=8,
                        symbol='circle'),
                    )
    )
    fig.update_layout(
        xaxis = dict(
            range = [0, t[-1]],
            title = {'text': 'Proper-Time/M'}
        ),
        yaxis = dict(
            range = [1.05*min(H), 1.05*max(H)],
            title = {'text': Htype}
        ),
        margin=dict(l=20, r=20, t=20, b=20),
        showlegend=False,        
    )
    return fig

def create_gw_figures(t, r, phi):
    Iddot = zwoc.get_Iddot(t, r, phi)
    # the derivatives make the Iddot smaller... pad them out with zeros
    Hplus = np.concatenate( (Iddot[0][0], [0,0,0,0]) )
    Hcross = np.concatenate( (Iddot[1][0], [0,0,0,0]) )
    # create the figures
    gw_plus_fig = makefig_gw(t, Hplus, "H+")
    gw_cross_fig = makefig_gw(t, Hcross, "Hx")
    return gw_plus_fig, gw_cross_fig, t, Hplus, Hcross

def create_effective_potential_figure(ell, E):
    r = np.linspace(rmin, rmax, Npoints_r)
    # when using the orbitcalc module, first set pars
    zwoc.G = G
    zwoc.M = M
    zwoc.ell = ell
    zwoc.E = E
    # calculate the potential
    V = zwoc.Veff(r)
    # calculate the Newtonian potential (not used because plot didn't work yet)
    VN = zwoc.VNeff(r)
    # find the positions two circular orbits (inner-unstable, outer-stable)
    r_inner, r_outer = zwoc.get_Veff_maxmin_r_values()
    Vmin = zwoc.Veff(r_outer)
    Vmax = zwoc.Veff(r_inner)
    deltaV = Vmax-Vmin
    # get info for nice plotting regions
    if (ell < 4*G*M):
        # find other radius where Veff = Vmax
        rwell = zwoc.get_rwell(Vmax)
    else:
        # find zero-crossing point of Veff
        zerocross = zwoc.get_rwell(Vmax)
        rwellA = 0.95*zerocross
        # set other one to rwellA plus 3 times dist to Vmin
        rwellB = 10*(r_outer-zerocross) + rwellA
    #
    #--- check the energy value to make sure it is valid
    #
    if E < Vmin:
        # just set to Vmin (stable circular orbit)
        E = Vmin
        zwoc.E = E
    elif E > Vmax:
        # forbid plunging orbits and hyperbolic orbits...
        # set to just below peak or just below zero
        #
        # FIXME: allow hyperbolic/plunging?
        #
        E1 = energy_frac_of_Vmax*deltaV + Vmin
        E2 = energy_frac_of_Vmax*(0 - Vmin) + Vmin
        E = min(E1, E2)
        zwoc.E = E
    #
    #--- Get peri and apoapsis
    #
    rp, ra = zwoc.get_peri_and_apoapsis()
    #
    #--- set the plotted window
    #
    #plot_rmax = 2*r_outer # rmax is twice circular orbit?
    if (zwoc.ell < 4*G*M):
        plot_rmin = 0.8*r_inner
        plot_rmax =1.02*rwell
        plot_Vmin = Vmin - deltaV*0.2
        plot_Vmax = Vmax + deltaV*0.2
    else:
        plot_rmin = 0.98*rwellA
        plot_rmax = 1.02*rwellB
        plot_Vmin = -1.2*np.abs(Vmin)
        plot_Vmax = 0.2*np.abs(Vmin)        
        #plot_rmin = 0
        #plot_rmax = 1.02*ra # rmax is apoapsis?
    #
    #--- create the figure
    #
    effpot_fig = makefig_effpot(
        r, V, VN, E, plot_rmin, plot_rmax, plot_Vmin, plot_Vmax
        )
    return effpot_fig, E
    
###############################
#    end helper functions     #
###############################

###############################
#         callbacks           #
###############################

#--- callback on page load 
@app.callback(
    Output('page-content', 'children'),
    Input('url', 'pathname')
)
def display_page(pathname):
    if (pathname == "/"):
        return initial_dashboard_page

#--- callback to open modal info window
@app.callback(
    Output("modal", "is_open"),
    [
        Input("open", "n_clicks"),
        Input("close", "n_clicks")
    ],
    [
        State("modal", "is_open")
    ],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

#--- client-side callback to "evolve" data
# javascript function within the clientside callback comment
app.clientside_callback(
    """
    function (n_intervals, data, offset, energy, gwdata) {
       if ( data == 0 ){
            return window.dash_clientside.no_update
       }
       offset = offset % data.r.length;
       const end = Math.min((offset + 1), data.r.length);
       const Earr = new Array(energy, energy);
       // note that we indicate the second trace "[1]" should be updated
       // and one point "1" is updated
       return [
          [
            {
               r: [data.r.slice(offset, offset+1)], 
               theta: [data.phi.slice(offset, offset+1)]
            },
            [1], 1
          ], 
          [
            {
               x: [data.r.slice(offset, offset+1)], 
               y: [Earr]
            },
            [2], 1
          ], 
          [
            {
               x: [gwdata.t.slice(offset, offset+1)], 
               y: [gwdata.plus.slice(offset, offset+1)],
            },
            [1], 1
          ], 
          [
            {
               x: [gwdata.t.slice(offset, offset+1)], 
               y: [gwdata.cross.slice(offset, offset+1)],
            },
            [1], 1
          ], 
          end
       ]
    }
    """, 
    [
        Output('orbitgraph', 'extendData'),
        Output('potentialgraph', 'extendData'),
        Output('gw_plus_graph', 'extendData'),
        Output('gw_cross_graph', 'extendData'),                
        Output('stored-offset', 'data'),
    ],
    [
        Input('orbit-interval', 'n_intervals'),
    ],
    [
        State('stored-orbit', 'data'),
        State('stored-offset', 'data'),
        State('stored-energy', 'data'),
        State('stored-gw-data', 'data'),        
    ]
)
#
#--- callback to re-make the effective potential plot
#
#     (triggered by any change to parameters (E, ell, e, rp),
#      or pressing the default button
#
@app.callback(
    [
        Output('potentialgraph', 'figure'),
        Output('angmom-val-str', 'value'),
        Output('energy-val-str', 'value'),
        Output('ecc-val-str', 'value'),
        Output('periap-val-str', 'value'),
        Output('stored-angmom', 'data'),
        Output('stored-energy', 'data'),        
        Output('stored-ecc', 'data'),
        Output('stored-periap', 'data'),        
    ],
    [
        Input('angmom-val-str', 'value'),
        Input('energy-val-str', 'value'),
        Input('ecc-val-str', 'value'),
        Input('periap-val-str', 'value'),
        Input('default-button', 'n_clicks'),
    ],
    [
        State('stored-angmom', 'data'),
        State('stored-energy', 'data'),
        State('stored-ecc', 'data'),
        State('stored-periap', 'data'),                
    ],
)
def remake_effective_potential(angmom_str, energy_str, ecc_str, periap_str,
                               n_clicks, ell_old, E_old, ecc_old, periap_old):
    # convert input strings to numbers (and check for invalid/erroneous)
    [ell_new, E_new, ecc_new, periap_new] = \
        get_all_values_from_strings(angmom_str, energy_str, ecc_str, periap_str, checkvalid=True)
    #print("--- input values:")
    #print("ell =", ell_new)
    #print("E =", E_new)          
    #print("ecc =", ecc_new)
    #print("periap =", periap_new)          
    # get trigger to see which change was made
    trigger = dash.callback_context.triggered[0]
    if ( (None in [ell_new, E_new, ecc_new, periap_new]) ):
        #
        # invalid/erroneous input -> revert to current values
        #
        ell = ell_old
        E = E_old
        ecc = ecc_old
        periap = periap_old
    elif ('default-button' in trigger['prop_id']):
        #
        # Default button pushed -> revert to default values
        #
        [ell, E, ecc, periap] = \
            get_all_values_from_strings(default_angmom_str, default_energy_str,
                                        default_ecc_str, default_periap_str)
    elif ( ('angmom' in trigger['prop_id']) | ('energy' in trigger['prop_id'])
           | (trigger['prop_id'] == '.') ):
           #
           # ang-mom/energy changed (or first run)
           #
           ell = ell_new
           E = E_new
           # get eccentricity and periapsis, adjusting energy if necessary
           [ecc, periap, E] = get_ecc_periap_and_check_E(ell, E)
    elif ('ecc' in  trigger['prop_id']):
        ecc = ecc_new
        periap = periap_new
        [ell, E, ecc, periap] = get_ell_E_and_check_pars(ecc, periap, version='ecc')
        if (None in [ell, E, ecc, periap]):
            ell = ell_old
            E = E_old
            ecc = ecc_old
            periap = periap_old
    elif ('periap' in  trigger['prop_id']):
        ecc = ecc_new
        periap = periap_new
        [ell, E, ecc, periap] = get_ell_E_and_check_pars(ecc, periap, version='periap')
        if (None in [ell, E, ecc, periap]):
            ell = ell_old
            E = E_old
            ecc = ecc_old
            periap = periap_old
    #print("--- output values:")
    #print("ell =", ell)
    #print("E =", E)          
    #print("ecc =", ecc)
    #print("periap =", periap)
    #
    #--- Create the effective potential figure and proceed...
    #
    effpot_fig, E = create_effective_potential_figure(ell, E)
    return [effpot_fig, str(ell), str(E), str(ecc), str(periap),
            ell, E, ecc, periap]

#
#--- callback to re-calculate orbit, save new data set, and restart figures
#
#    (triggered by stored-energy from "remake_effective_potential" callback)
#
@app.callback(
    [
        Output('orbitgraph', 'figure'),
        Output('stored-orbit', 'data'),
    ],
    [
        Input('stored-energy', 'data')
    ],
    [
        State('angmom-val-str', 'value'),        
        State('energy-val-str', 'value'),
        State('ecc-val-str', 'value'),
        State('periap-val-str', 'value'),        
    ],
)
def recalculate_orbit(energy, angmom_str, energy_str, ecc_str, periap_str):
    # The input energy and angular momentum values should
    # always be valid, because if a user entry is not valid,
    # the remake_effective_potential() callback will adjust
    # them to their (valid) default values.    
    [ell, E, ecc, periap] = \
        get_all_values_from_strings(angmom_str, energy_str, ecc_str, periap_str)
    fig_orb, t, r_t, phi_t = create_orbit_figure(ell, E)
    stored_data = dict(t = t, r = r_t, phi = phi_t,
                       resolution=default_orbit_resolution)
    return fig_orb, stored_data

#
#--- callback to re-calculate gravitational wave data on orbit update
#
#    (triggered by stored-orbit from "recalculate_orbit" callback)
#
@app.callback(
    [
        Output('gw_plus_graph', 'figure'),
        Output('gw_cross_graph', 'figure'),
        Output('stored-gw-data', 'data'),        
    ],
    [
        Input('stored-orbit', 'data'),
    ],
    [
        State('stored-energy', 'data'),
        State('stored-angmom', 'data'),        
        State('stored-ecc', 'data'),
        State('stored-periap', 'data'),        
    ],
)
def recalculate_gw(orbit_data, E, ell, ecc, periap):
    t = orbit_data['t']
    r = orbit_data['r']
    phi = orbit_data['phi']
    gw_plus_fig, gw_cross_fig, t, Hp, Hc =  create_gw_figures(t, r, phi)
    stored_data = dict(t = t, plus = Hp, cross = Hc,
                       resolution=default_orbit_resolution)
    # output the data to user (this should usually be turned off!)
    if print_out_data:
        print("################################")
        print("# orbit and GW data for:")
        print("#    E =", E, "l =", ell)
        print("################################")
        print("# t  r(t)  phi(t)  Hplus(t)  Hcross(t)")
        for i in range(len(t)):
            print(t[i], r_t[i], phi_t[i], Hp[i], Hc[i])    
    return gw_plus_fig, gw_cross_fig, stored_data

###############################
#        end callbacks        #
###############################

###############################
#        main routine         #
###############################

#
#--- Create initial orbit figures
#
ell = get_number_from_string(default_angmom_str)
E = get_number_from_string(default_energy_str)
init_orbit_fig, t, r_t, phi_t = create_orbit_figure(ell, E)
init_orbit_data = dict(r = r_t, phi = phi_t,
                       resolution = default_orbit_resolution)
init_pot_fig, E = create_effective_potential_figure(ell, E)
init_gw_plus_fig, init_gw_cross_fig, t, Hp, Hc = \
    create_gw_figures(t, r_t, phi_t)
init_gw_data = dict(t = t, plus = Hp, cross = Hc,
                    resolution=default_orbit_resolution)                    
#
#--- Grab the webpage from the zoom-whirl app html module
#
initial_dashboard_page = zwah.make_dashboard_webpage(
    #None, None, None,
    init_orbit_fig, init_orbit_data, init_pot_fig,
    init_gw_cross_fig, init_gw_plus_fig, init_gw_data,
    default_angmom_str,
    default_energy_str,
    default_ecc_str,
    default_periap_str,    
    )

# run the app
if __name__ == '__main__':
    app.run_server(debug=True)
    
