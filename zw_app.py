import dash
import dash.html as html
import dash.dcc as dcc
import dash_bootstrap_components as dbc
import dash.exceptions as dex
import plotly.express as px
import plotly.graph_objects as pgo
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_pdf import FigureCanvasPdf
import pandas as pd
import datetime as dt
import zipfile
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
default_tmax_str = str(zwoc.tf_default) # 1000
default_speed_str = "4" #

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

def create_orbit_figure(ell, E, tmax):
    # evolve the orbit
    t, r_t, phi_t, rp, ra, ecc = zwoc.get_orbit(ell, E, tmax)
    # make a figure
    fig = makefig_orbit(r_t, phi_t, ra)
    # return the figure and the data
    return fig, t, r_t, phi_t

def ell_of_er(e, rp):
    """ 
    mapping (e, rp) -> (ell, E)
    """
    return (1 + e) * rp / np.sqrt( (1 + e) * rp - (3 + e**2) )

def E_of_er(e, rp):
    """ 
    mapping (e, rp) -> (ell, E)
    """
    return ( (1 - e) / (2 * rp) * (4 - (1 + e) * rp)
             / ( (1 + e) * rp - (3 + e**2) ) )
        
def get_ell_E_and_check_pars(ecc, periap, version=None):
    """
    Convert (e, rp) -> (ell, E), and verify that the result is valid.
    If the result is not valid, try to change either eccentricity or
    periapsis slightly (whichever was not set by user) to make it so.
    
    The "version" specifies which parameter was altered by the user. 
    If eccentricity was altered, then it is assumed fixed and periapsis
    is adjusted, if necessary, to get a valid output.  And vice versa.
    """
    if version == 'ecc':
        rp = periap
        #
        # Check validity:
        #
        #   It turns out that the condition for the ell(e, rp)
        #   and E(e, rp) equations to be true, i.e., that r0 < rp
        #   (see Poisson & Will, Sec 5.6.3):
        #
        #          rp > 2*(3+e)/(1+e)
        #
        #   is sufficient (see my 2023-02-18 notes) to make both
        #   the denominator of [ell(e, rp)]^2 valid, i.e.,
        #
        #         (1+e)*rp - (3+e^2) > 0
        #
        #   and for the condition that there is in fact a
        #   "potential well":
        #
        #          ell^2 > 12
        #
        #   So, we need only check that top requirement, and, if
        #   it is not satisfied, change either the periapsis
        #   or eccentricity (whichever was *not* set by user)
        #   to fit it.
        #
        #   But our additional maximum value of ell places one
        #   more non-physical requirement on (e, rp).
        #
        # check requirement and adjust periapsis if necessary:
        if (rp < 2*(3+ecc)/(1+ecc)) :
            rp = 1.05 * (2*(3+ecc)/(1+ecc))
        # set values of ell(e,rp) and E(e,rp)
        ell = ell_of_er(ecc, rp)
        E = E_of_er(ecc, rp)
        # verify that angular momentum is not too large
        if ( ell > ell_max ):
            # failure: will revert to prior values
            rp = None
        return [ell, E, ecc, rp]
    elif version == 'periap':
        e = ecc
        # check requirement and adjust eccentricity if necessary:
        if (periap < 2*(3+e)/(1+e)) :
            e = np.max([0, 1.05* (6 - periap)/(periap - 2)])
        # set values of ell(e,rp) and E(e,rp)
        ell = ell_of_er(e, periap)
        E = E_of_er(e, periap)
        # verify that angular momentum is not too large
        if ( ell > ell_max ):
            # failure: will revert to prior values
            e = None
        return [ell, E, e, periap]

def get_ecc_periap_and_check_E(ell, E):
    # find the positions two circular orbits (inner-unstable, outer-stable)
    r_inner, r_outer = zwoc.get_Veff_maxmin_r_values(ell)
    Vmin = zwoc.Veff(r_outer, ell)
    Vmax = zwoc.Veff(r_inner, ell)
    deltaV = Vmax-Vmin
    #
    #--- check the energy value to make sure it is valid
    #
    if E < Vmin:
        # just set to Vmin (stable circular orbit)
        E = Vmin
    elif E > Vmax:
        # forbid plunging orbits and hyperbolic orbits...
        # set to just below peak or just below zero
        #
        # FIXME: allow hyperbolic/plunging?
        #
        E1 = energy_frac_of_Vmax*deltaV + Vmin
        E2 = energy_frac_of_Vmax*(0 - Vmin) + Vmin
        E = min(E1, E2)
    # get periapsis and apoapsis and calculate eccentricity
    rp, ra = zwoc.get_peri_and_apoapsis(ell, E)
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
                     title = {'text': 'radius, r (M)'}),
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
            title = {'text': 'proper time, \u03c4 (M)'}
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
    # calculate the potential
    V = zwoc.Veff(r, ell)
    # calculate the Newtonian potential (not used because plot didn't work yet)
    VN = zwoc.VNeff(r, ell)
    # find the positions two circular orbits (inner-unstable, outer-stable)
    r_inner, r_outer = zwoc.get_Veff_maxmin_r_values(ell)
    Vmin = zwoc.Veff(r_outer, ell)
    Vmax = zwoc.Veff(r_inner, ell)
    deltaV = Vmax-Vmin
    # get info for nice plotting regions
    if (ell < 4*G*M):
        # find other radius where Veff = Vmax
        rwell = zwoc.get_rwell(Vmax, ell)
    else:
        # find zero-crossing point of Veff
        zerocross = zwoc.get_rwell(Vmax, ell)
        rwellA = 0.95*zerocross
        # set other one to rwellA plus 3 times dist to Vmin
        rwellB = 10*(r_outer-zerocross) + rwellA
    #
    #--- check the energy value to make sure it is valid
    #
    if E < Vmin:
        # just set to Vmin (stable circular orbit)
        E = Vmin
    elif E > Vmax:
        # forbid plunging orbits and hyperbolic orbits...
        # set to just below peak or just below zero
        #
        # FIXME: allow hyperbolic/plunging?
        #
        E1 = energy_frac_of_Vmax*deltaV + Vmin
        E2 = energy_frac_of_Vmax*(0 - Vmin) + Vmin
        E = min(E1, E2)
    #
    #--- Get peri and apoapsis
    #
    rp, ra = zwoc.get_peri_and_apoapsis(ell, E)
    #
    #--- set the plotted window
    #
    #plot_rmax = 2*r_outer # rmax is twice circular orbit?
    if (ell < 4*G*M):
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
    rind = np.where((r > plot_rmin) & (r < plot_rmax))
    E_v_r = E*(r*0 + 1)
    return effpot_fig, E, r[rind], V[rind], E_v_r[rind]

def create_plot_for_export(bytes_io, ep_r, ep_V, ep_E, orb_phi, orb_r,
                           gw_t, gw_Hplus, gw_Hcross,
                           ell, E, ecc, periap, M_over_m):
    """
    This function creates a matplotlib plot *without* using pyplot, which 
    generates this error when used:
      "Starting a Matplotlib GUI outside of the main thread will likely fail."
    I found the solution at:
            https://stackoverflow.com/questions/69650149    
    """
    # use gridspec to make 4 figures in nice aspect ratio
    gs = mpl.gridspec.GridSpec(7, 2) # nrows, ncols
    # call the "Figure" method directly from matplotlib
    fig = mpl.figure.Figure(figsize=(8, 10.5))
    # effective potential plot
    ax1 = fig.add_subplot(gs[0:3,0])
    ax1.plot(ep_r, ep_V)
    ax1.plot(ep_r, ep_E)
    ax1.set_xlabel(r'$r/M$')
    ax1.set_ylabel(r'$V_{\rm eff}$')
    ax1.text(0.8, 0.2, r"$\ell$" + f" = {ell:.6f}\n" + r"$E$" + f" = {E:.6f}",
             horizontalalignment='center',
             verticalalignment='center', transform=ax1.transAxes)
    # orbit plot
    ax3 = fig.add_subplot(gs[0:3,1], polar=True)
    ax3.plot(orb_phi, orb_r)
    ax3.text(0.8, 0, r"$e$" + f" = {ecc:.3f}\n" + r"$r_p$" + f" = {periap:.2f}",
             horizontalalignment='center',
             verticalalignment='top', transform=ax3.transAxes)
    # Gwave-Hplus plot
    ax2 = fig.add_subplot(gs[3:5,:])
    ax2.plot(gw_t, gw_Hplus)
    ax2.set_xlabel(r'$\tau/M$')
    ax2.set_ylabel(r'$H_+$')
    ax2.set_title(f"M/m = {round(M_over_m):d}")
    # Gwave-Hcross plot
    ax4 = fig.add_subplot(gs[5:7,:])
    ax4.plot(gw_t, gw_Hcross)
    ax4.set_xlabel(r'$\tau/M$')
    ax4.set_ylabel(r'$H_x$')
    # gather plots and adjust margins
    fig.tight_layout()
    # write the pdf image to the buffer that came in as an argument
    #  (again, without using pyplot's "savefig")
    canvas = FigureCanvasPdf(fig)
    canvas.print_pdf(bytes_io)

def create_dataframes_of_stored_data(orbit_data, gw_data, effpot_data,
                                 ell, E, ecc, periap):
    # create dataframes of stored data
    df_orb = pd.DataFrame({
        't': orbit_data['t'],
        'r': orbit_data['r'],
        'phi': orbit_data['phi']})
    df_gw = pd.DataFrame({
        't': gw_data['t'],
        'Hplus': gw_data['plus'],
        'Hcross': gw_data['cross']})
    df_effpot = pd.DataFrame({
        'r': effpot_data['r'],
        'Veff': effpot_data['Veff'],
        'E': effpot_data['E'],})
    # create dataframe of parameter values
    df_pars = pd.DataFrame({
        'ell': ell,
        'E': E,
        'ecc': ecc,
        'periap': periap,
        'M_over_m': [zwoc.M/zwoc.m],})
    return [df_orb, df_gw, df_effpot, df_pars]


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
    function (n_intervals, data, energy, gwdata, offset, offset_speed) {
       if ( data == 0 ){
            return window.dash_clientside.no_update
       }
       offset = offset % data.r.length;
       const end = Math.min((offset + offset_speed), data.r.length);
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
        State('stored-energy', 'data'),
        State('stored-gw-data', 'data'),        
        State('stored-offset', 'data'),
        State('stored-speed', 'data'),        
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
        Output('tmax-val-str', 'value'),
        Output('stored-angmom', 'data'),
        Output('stored-energy', 'data'),        
        Output('stored-ecc', 'data'),
        Output('stored-periap', 'data'),
        Output('stored-tmax', 'data'),        
        Output('stored-effpot', 'data'),
    ],
    [
        Input('angmom-val-str', 'value'),
        Input('energy-val-str', 'value'),
        Input('ecc-val-str', 'value'),
        Input('periap-val-str', 'value'),
        Input('tmax-val-str', 'value'),
        Input('default-button', 'n_clicks'),
    ],
    [
        State('stored-angmom', 'data'),
        State('stored-energy', 'data'),
        State('stored-ecc', 'data'),
        State('stored-periap', 'data'),
        State('stored-tmax', 'data'),
        State('stored-offset', 'data'),                
    ],
)
def remake_effective_potential(angmom_str, energy_str, ecc_str, periap_str, tmax_str,
                               n_clicks, ell_cur, E_cur, ecc_cur, periap_cur, tmax_cur, offset):
    # convert input strings to numbers (and check for invalid/erroneous)
    [ell_new, E_new, ecc_new, periap_new] = \
        get_all_values_from_strings(angmom_str, energy_str, ecc_str, periap_str, checkvalid=True)
    # get trigger to see which change was made
    trigger = dash.callback_context.triggered[0]
    # assume tmax unchanged for now
    tmax = int(tmax_cur)
    if ( (None in [ell_new, E_new, ecc_new, periap_new]) ):
        #
        # invalid/erroneous input -> revert to current values
        #
        ell = ell_cur
        E = E_cur
        ecc = ecc_cur
        periap = periap_cur
    elif ('default-button' in trigger['prop_id']):
        #
        # Default button pushed -> revert to default values
        #
        [ell, E, ecc, periap] = \
            get_all_values_from_strings(default_angmom_str, default_energy_str,
                                        default_ecc_str, default_periap_str)
        # and reset tmax to default value
        tmax = int(zwoc.tf_default)
    elif ( ('angmom' in trigger['prop_id']) | ('energy' in trigger['prop_id'])
           | (trigger['prop_id'] == '.') ):
           #
           # ang-mom/energy changed (or first run)
           #
           ell = ell_new
           E = E_new
           # get eccentricity and periapsis, adjusting energy if necessary
           [ecc, periap, E] = get_ecc_periap_and_check_E(ell, E)
    elif ('ecc' in trigger['prop_id']):
        ecc = ecc_new
        periap = periap_new
        [ell, E, ecc, periap] = get_ell_E_and_check_pars(ecc, periap, version='ecc')
        if (None in [ell, E, ecc, periap]):
            ell = ell_cur
            E = E_cur
            ecc = ecc_cur
            periap = periap_cur
    elif ('periap' in trigger['prop_id']):
        ecc = ecc_new
        periap = periap_new
        [ell, E, ecc, periap] = get_ell_E_and_check_pars(ecc, periap, version='periap')
        if (None in [ell, E, ecc, periap]):
            ell = ell_cur
            E = E_cur
            ecc = ecc_cur
            periap = periap_cur
    elif ('tmax-val' in trigger['prop_id']):
        # all parameters are unchanged
        ell = ell_cur
        E = E_cur
        ecc = ecc_cur
        periap = periap_cur
        # but change the time of integration
        tmax = int(get_number_from_string(tmax_str))
        # as long as it is within the allowed range
        if ( (tmax >= float(zwah.tmax_min_str))
             & (tmax <= float(zwah.tmax_max_str)) ):
            # ok, will store new tmax
            pass
        else:
            # if out of bounds, leave tmax unchanged
            tmax = int(tmax_cur)
    #
    #--- Create the effective potential figure and proceed...
    #
    effpot_fig, E, r, V, E_v_r = create_effective_potential_figure(ell, E)
    stored_data = dict(r = r, Veff = V, E = E_v_r)
    return [effpot_fig, str(ell), str(E), str(ecc), str(periap), str(tmax),
            ell, E, ecc, periap, tmax, stored_data]

#
#--- callback to change speed of orbit
#
#     (triggered by changing the speed (in "x"))
#
@app.callback(
    [
        Output('stored-speed', 'data'),
        Output('speed-val-str', 'value')
    ],
    [
        Input('speed-val-str', 'value'),
    ],
    [
        State('stored-speed', 'data')        
    ],
)
def change_speed(speed_str, oldspeed):
    newspeed = int(get_number_from_string(speed_str))
    if ( (newspeed >= float(zwah.speed_min_str))
         & (newspeed <= float(zwah.speed_max_str)) ):
        return newspeed, str(newspeed)
    else:
        return oldspeed, str(oldspeed)

#
#--- callback to download a nice matplotlib pdf of current plots
#
#     (triggered by pressing the "export plot" button)
#
@app.callback(
    [
        Output('plot-download', 'data')
    ],
    [
        Input('plot-button', 'n_clicks'),
    ],
    [
        State('stored-orbit', 'data'),
        State('stored-gw-data', 'data'),
        State('stored-effpot', 'data'),            
        State('stored-angmom', 'data'),
        State('stored-energy', 'data'),
        State('stored-ecc', 'data'),
        State('stored-periap', 'data'),                
    ],
)
def download_plot(n_clicks, orbit_data, gw_data, effpot_data,
                  ell, E, ecc, periap):
    # create (dated and parameter-labeled) filename for plot pdf
    nowstr = dt.datetime.now().strftime("%Y-%m-%d_%H%M-%S")
    plot_filename = "zoomwhirl-data_L_" + f"{ell:.3e}" + "_E_" \
        + f"{E:.3e}" +    "_" + nowstr + ".pdf"
    # this function sends the destination "bytes_io" to a
    # subroutine that creates a plot with matplotlib,
    # writing the pdf output to "bytes_io"
    def write_plot(bytes_io):
        create_plot_for_export(bytes_io,
            effpot_data['r'], effpot_data['Veff'], effpot_data['E'],
            orbit_data['phi'], orbit_data['r'],
            gw_data['t'], gw_data['plus'], gw_data['cross'],
            ell, E, ecc, periap, zwoc.M/zwoc.m)
    return [dcc.send_bytes(write_plot, plot_filename)]

#
#--- callback to download all data sets currently plotted
#
#     (triggered by pressing the "download data" button)
#
@app.callback(
    [
        Output('data-download', 'data')
    ],
    [
        Input('download-button', 'n_clicks'),
    ],
    [
        State('stored-orbit', 'data'),
        State('stored-gw-data', 'data'),
        State('stored-effpot', 'data'),            
        State('stored-angmom', 'data'),
        State('stored-energy', 'data'),
        State('stored-ecc', 'data'),
        State('stored-periap', 'data'),                
    ],
)
def download_data(n_clicks, orbit_data, gw_data, effpot_data, ell, E, ecc, periap):
    # put the stored data into dataframes
    [df_orb, df_gw, df_effpot, df_pars] = \
        create_dataframes_of_stored_data(orbit_data, gw_data, effpot_data,
                                         ell, E, ecc, periap)
    # convert all dataframes into csv-file strings for output to zip
    dfs = [df_orb, df_gw, df_effpot, df_pars]
    df_names = ['orbit', 'grav-wave', 'eff-pot', 'pars']
    df_strs = []
    for df in dfs:
        dfstr = df.to_csv(index=False)
        df_strs.append(dfstr)
    # create (dated and parameter-labeled) filename for zip file
    nowstr = dt.datetime.now().strftime("%Y-%m-%d_%H%M-%S")
    zip_filename = "zoomwhirl-data_L_" + f"{ell:.3e}" + "_E_" \
        + f"{E:.3e}" +    "_" + nowstr + ".zip"
    # This function takes the bytes_io argument as the place
    # to write the zip file.
    #    (found here: https://stackoverflow.com/questions/67917360)
    def write_archive(bytes_io):
        with zipfile.ZipFile(bytes_io, mode="w") as zf:
            # add all dataset csv files
            for df, n in zip(df_strs, df_names):
                zf.writestr(n + '.csv', df)
            # add in the python (matplotlib) script that creates pretty plots
            zf.write("additional-code/plot-output-data.py", arcname="plot-output-data.py")
            # add instructions regarding the files
            zf.write("additional-code/plot-output-data_README.txt", arcname="README")
    return [dcc.send_bytes(write_archive, zip_filename)]

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
        State('stored-tmax', 'data'),        
    ],
)
def recalculate_orbit(energy, angmom_str, energy_str, ecc_str, periap_str, tmax):
    # The input energy and angular momentum values should
    # always be valid, because if a user entry is not valid,
    # the remake_effective_potential() callback will adjust
    # them to their (valid) default values.    
    [ell, E, ecc, periap] = \
        get_all_values_from_strings(angmom_str, energy_str, ecc_str, periap_str)
    fig_orb, t, r_t, phi_t = create_orbit_figure(ell, E, tmax)
    stored_data = dict(t = t, r = r_t, phi = phi_t)
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
    stored_data = dict(t = t, plus = Hp, cross = Hc)
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
tmax = get_number_from_string(default_tmax_str)
init_orbit_fig, t, r_t, phi_t = create_orbit_figure(ell, E, tmax)
init_orbit_data = dict(r = r_t, phi = phi_t)
init_pot_fig, E, r, V, E_v_r = create_effective_potential_figure(ell, E)
init_gw_plus_fig, init_gw_cross_fig, t, Hp, Hc = \
    create_gw_figures(t, r_t, phi_t)
init_gw_data = dict(t = t, plus = Hp, cross = Hc)
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
    default_tmax_str,
    default_speed_str,        
    )

# run the app
if __name__ == '__main__':
    app.run_server(debug=True)
    
