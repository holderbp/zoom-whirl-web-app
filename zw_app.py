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
#
#--- effective potential graph
#
# Npoints, and min and max values for making graph
Npoints_r = 1000
rmin = 0.1
rmax = 100
# max allowed ang momentum value
ell_max = 10.0
# for choosing an energy just below Vmax
energy_frac_of_Vmax = 0.9
# starting values
default_angmom_str = "3.7" # "4.0"
default_energy_str = "-0.0355" # "-0.03"

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
        margin=dict(l=20, r=20, t=20, b=20),
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

def get_number_from_string(valstr):
    try:
        val = float(valstr)
        return val
    except ValueError:
        return None

def angmom_is_valid(L):
    # Disallow cases where there are no bound orbits
    if ( (L**2 > 12*(G*M)**2) & (L < ell_max) ):
        return True
    else:
        return False

def makefig_effpot(r, V, E, win_rmin, win_rmax, win_Vmin, win_Vmax):
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
                     title = {'text': 'Radius'}),
        yaxis = dict(range = [win_Vmin, win_Vmax],
                     title = {'text': 'Veff'}),
        margin=dict(l=20, r=20, t=20, b=80),
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
            title = {'text': 'time'}
        ),
        yaxis = dict(
            title = {'text': Htype}
        ),
        margin=dict(l=20, r=20, t=20, b=80),
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
    #
    #--- Get peri and apoapsis
    #
    rp, ra = zwoc.get_peri_and_apoapsis()
    #
    #--- set the plotted window
    #
    plot_rmin = 0
    #plot_rmax = 2*r_outer # rmax is twice circular orbit?
    plot_rmax = 1.02*ra # rmax is apoapsis?
    plot_Vmin = Vmin - deltaV*0.2
    plot_Vmax = Vmax + deltaV*0.2
    #
    #--- create the figure
    #
    effpot_fig = makefig_effpot(
        r, V, E, plot_rmin, plot_rmax, plot_Vmin, plot_Vmax
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
#     (whenever angular momentum and E are adjusted)
#
@app.callback(
    [
        Output('potentialgraph', 'figure'),
        Output('angmom-val-str', 'value'),
        Output('energy-val-str', 'value'),
        Output('stored-energy', 'data'),
    ],
    [
        Input('angmom-val-str', 'value'),
        Input('energy-val-str', 'value'),
    ],
    [
    ],
)
def remake_effective_potential(angmom_str, energy_str):
    ell = get_number_from_string(angmom_str)
    E = get_number_from_string(energy_str)
    #
    #--- If either of these strings are not valid numbers
    #    revert to default values and make original figure
    #
    #--- Do the same if the angular momentum falls outside
    #    the allowed range
    #
    if ( (ell is None) | (E is None)
         | (not angmom_is_valid(ell)) ):
        ell = get_number_from_string(default_angmom_str)
        E = get_number_from_string(default_energy_str)        
        #
        # FIXME:
        #
        #   * should probably have an Alert signal to user here,
        #     explaining why their update failed
        #
    #
    #--- otherwise calculate the new effective potential using ell
    #      (and get new energy if found to be invalid)
    #
    effpot_fig, E = create_effective_potential_figure(ell, E)
    return [effpot_fig, str(ell), str(E), E]

#--- callback to re-calculate orbit, save new data set, and restart figures
@app.callback(
    [
        Output('orbitgraph', 'figure'),
        Output('stored-orbit', 'data'),
    ],
    [
        Input('recalculate-button', 'n_clicks'),
        Input('stored-energy', 'data')
    ],
    [
        State('energy-val-str', 'value'),
        State('angmom-val-str', 'value'),        
    ],
)
def recalculate_orbit(n_clicks, energy, energy_str, angmom_str):
    # The input energy and angular momentum values should
    # always be valid, because if a user entry is not valid,
    # the remake_effective_potential() callback will adjust
    # them to their (valid) default values.
    ell = get_number_from_string(angmom_str)
    E = get_number_from_string(energy_str)
    fig_orb, t, r_t, phi_t = create_orbit_figure(ell, E)
    stored_data = dict(t = t, r = r_t, phi = phi_t,
                       resolution=default_orbit_resolution)
    return fig_orb, stored_data

#--- callback to re-calculate gravitational wave data on orbit update
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
    ],
)
def recalculate_gw(orbit_data):
    t = orbit_data['t']
    r = orbit_data['r']
    phi = orbit_data['phi']
    gw_plus_fig, gw_cross_fig, t, Hp, Hc =  create_gw_figures(t, r, phi)
    stored_data = dict(t = t, plus = Hp, cross = Hc,
                       resolution=default_orbit_resolution)
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
    default_energy_str
    )

# run the app
if __name__ == '__main__':
    app.run_server(debug=True)
    
