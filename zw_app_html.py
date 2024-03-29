import dash.html as html
import dash.dcc as dcc
import dash_bootstrap_components as dbc
import zw_app_orbitcalc as zwoc

effpot_fontsize = '10pt'
effpot_fontsize_small = '8pt'

height_header = '5vh'
height_firstrow = '50vh'
height_secondrow = '35vh'
height_secondrow_card = '45vh'

# max allowed angular momentum 
ell_max_str = "20" # M (converted to float and used by the app)
tmax_min_str = "100"
tmax_max_str = "10000"
speed_min_str = "1"
speed_max_str = "10"
bhmass_min_str = "10000"

###############################
#   dashboard webpage setup   #
###############################
def make_dashboard_webpage(
        init_orbit_fig, init_orbit_data,
        init_pot_fig,
        init_gw_cross_fig, init_gw_plus_fig,
        init_gw_data,
        default_angmom_str,
        default_energy_str,
        default_ecc_str,
        default_periap_str,
        default_tmax_str,
        default_speed_str,
        default_bhmass_str,         
):
    dashboard_page = html.Div(
        children=[
            #
            # =======================================
            # These items apply to all rows & columns
            # =======================================
            #
            # Downloads:
            #
            #    for downloading zipped file of all csv data frames
            #
            dcc.Download(id='data-download'),
            dcc.Download(id='plot-download'),            
            #
            # Intervals:
            #
            #    orbit-interval: for continuous clientside update (animation)
            # 
            dcc.Interval(id="orbit-interval",
                         interval=20, # microseconds
                         disabled=False,
                         n_intervals=0,
                         max_intervals=-1),
            #
            # Stores: client-stored json data sets:
            #
            #     offset
            #     orbit plot data
            #     energy value
            #     ang momentum value
            #     grav wave plot data
            #
            dcc.Store(id='stored-offset',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-angmom',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-energy',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-ecc',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-periap',
                      storage_type='memory',
                      data=0,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-orbit',
                      storage_type='memory',
                      data=0,
                      #data=init_orbit_data,                      
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-gw-data',
                      storage_type='memory',
                      data=0,
                      #data=init_gw_data,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-effpot',
                      storage_type='memory',
                      data=0,
                      #data=init_gw_data,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-tmax',
                      storage_type='memory',
                      data=int(zwoc.tf_default),
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-speed',
                      storage_type='memory',
                      data=1,
                      clear_data=False,
                      modified_timestamp=-1),
            dcc.Store(id='stored-bhmass',
                      storage_type='memory',
                      data=int(zwoc.M_over_m_default),
                      clear_data=False,
                      modified_timestamp=-1),
            #
            #=== Modal (information pop-up window)
            #
            # modal = html.Div(
            #     [
            #         dbc.Button("Open modal", id="open", n_clicks=0),
            #         dbc.Modal(
            #             [
            #                 dbc.ModalHeader(dbc.ModalTitle("Header")),
            #                 dbc.ModalBody("This is the content of the modal"),
            #                 dbc.ModalFooter(
            #                     dbc.Button(
            #                         "Close", id="close", className="ms-auto", n_clicks=0
            #                     )
            #                 ),
            #             ],
            #             id="modal",
            #             is_open=False,
            #         ),
            #     ]
            # )
            #
            #=== Header: Title, Information button, logos
            #
            html.Div( style = {'display': 'flex'}, children=[
                html.Div( style =
                          {
                              'width': '100%',
                              'margin-left': 25,
                              'margin-top': 10,
                              'margin-right': 10
                          },
                          children=[ 
                    html.Div( style={'float': 'left'}, children=[
                        html.H3("\"Zoom-Whirl\" orbits and their gravitational waves"),
                        html.H4("Exploring the geodesics of the Schwarzschild metric")
                        ]),
                    html.Div( style={'float': 'left',
                                     'margin-left': '5%',
                                     'font-size': effpot_fontsize,
                                     'textAlign': 'center'}, children=[
                        html.Div( children=[
                            html.Div( children="orbit speed: 1\u00d7 \u2264",  # \u2264 is ≤, \u00d7 is ×
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 5,
                                    'margin-bottom': 5,
                                    'font-size' : effpot_fontsize,
                                }
                            ),
                            dcc.Input(
                                id='speed-val-str', type='text',
                                value = default_speed_str,
                                placeholder = default_speed_str,
                                # don't allow update until user enters
                                debounce = True,
                                size = '2',
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 2,
                                    'font-size' : effpot_fontsize,
                                }
                            ),
                            html.Div(
                                children="\u00d7 \u2264 " + speed_max_str + "\u00d7",
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 5,
                                    'font-size' : effpot_fontsize,
                                }
                            ),
                        ]),
                        html.Div( children=[
                            html.Div(
                                children= "\u03c4\u2098\u2090\u2093 : " + tmax_min_str + " M \u2264",
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 5,
                                    'margin-bottom': 5,
                                    'font-size' : effpot_fontsize,
                                }
                            ),                                        
                            dcc.Input(
                                id='tmax-val-str', type='text',
                                value = default_tmax_str,
                                placeholder = default_tmax_str,
                                # don't allow update until user enters
                                debounce = True,
                                size = '7',
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 5,
                                    'font-size' : effpot_fontsize,
                                }
                            ),                                        
                            html.Div(
                                children="\u2264 " + tmax_max_str + " M",
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 5,
                                    'font-size' : effpot_fontsize,
                                }
                            ),
                        ]),
                        html.Div( children=[
                            html.Button('export data (zip)',
                                        id='download-button', n_clicks=0,
                                        style ={
                                            'font-size' : effpot_fontsize_small,
                                            'margin-bottom' : '5px',
                                        }),
                            html.Button('export plots (pdf)',
                                        id='plot-button', n_clicks=0,
                                        style ={
                                            'font-size' : effpot_fontsize_small,
                                        }),
                        ]),
                    ]),                        
                    html.Div( style={'float': 'left',
                                     'margin-top': '25px',
                                     'margin-left': '5%',
                                     'font-size': effpot_fontsize,
                                     'textAlign': 'center'}, children=[
                        html.A("arXiv", href="https://arxiv.org/abs/2303.04072",
                               target="_blank"),
                        html.Br(),
                        html.A("github", href="https://github.com/holderbp/zoom-whirl-web-app",
                               target="_blank"),
                    ]),
                    html.Div( style={'float': 'right'}, children=[
                        html.A(
                            html.Img(
                                src="assets/northwestern-ciera.png",
                                style = {'float': 'right', 'height': '35px', 'margin-top': '20px'}
                            ),
                            href="https://ciera.northwestern.edu",
                            target="_blank"),                                                       
                        html.A(
                            html.Img(
                                src="assets/gvsu.png",
                                style = {'float': 'right', 'height': '35px', 'margin-top': '20px',
                                         'margin-right': '20px'}
                            ),
                            href="https://www.gvsu.edu/physics/",
                            target="_blank"),                            
                        ]),
                    ]),
                ]),
            #
            #=== First Row
            #
            dbc.Row([
                dbc.CardGroup(style = {'height' : height_firstrow}, children=[
                    dbc.Card([
                        dbc.CardHeader("Orbit"),
                        dbc.Row(style={'height': '100%'}, children=[
                            dbc.Col(width=5, style={'margin-left': '20px'}, children=[
                                html.Div( children=[
                                    html.Br(),
                                    html.P("To generate a new orbit, specify:",
                                           style={
                                               'font-size': effpot_fontsize,
                                               'font-weight': 'bold',
                                               'color': '#636EFA',
                                               'margin-bottom': '5px',
                                           }),
                                    html.P("[Angular Momentum, Energy]",
                                           style={
                                               'margin-left': '25px',
                                               'font-size': effpot_fontsize,
                                               'margin-bottom': '5px',
                                               'margin-top': '5px',                                                                                          }),
                                    html.P("or:",
                                           style={
                                               'font-size': effpot_fontsize,
                                               'font-weight': 'bold',          
                                               'color': '#636EFA',
                                               'margin-bottom': '5px',
                                               'margin-top': '5px',
                                           }),
                                    html.P("[Eccentricity, Periapsis]",
                                           style={
                                               'margin-left': '25px',
                                               'font-size': effpot_fontsize,
                                               'margin-bottom': '5px',
                                               'margin-top': '5px',
                                           }),
                                    html.P("Orbits with energy near the peak of the effective potential are ZW-type.",
                                           style={
                                               'font-size': effpot_fontsize,
                                               'font-weight': 'bold',
                                               'color': '#636EFA',
                                               'margin-bottom': '10px',
                                               'margin-top': '5px',
                                           }),
                                    html.Div(style={
                                        'textAlign': 'center',
                                        'border': '2px lightgrey solid',
                                    }, children=[
                                        html.Label("Eccentricity: ",
                                                   style={
                                                       'font-size': effpot_fontsize,
                                                       'margin-top': '10px',
                                                       'margin-right': 10,
                                                   }),
                                        html.Div(
                                            children="0 \u2264",
                                            style={
                                                'display': 'inline-block',
                                                'margin-right': 5,
                                                'font-size' : effpot_fontsize,
                                            }
                                        ),                                        
                                        dcc.Input(
                                            id='ecc-val-str', type='text',
                                            value = default_ecc_str,
                                            placeholder = default_ecc_str,
                                            # don't allow update until user enters
                                            debounce = True,
                                            size = '7',
                                            style={
                                                'display': 'inline-block',
                                                'margin-right': 5,
                                                'font-size' : effpot_fontsize,
                                            }
                                        ),
                                        html.Div(
                                            children=" < 1",
                                            style={
                                                'display': 'inline-block',
                                                'font-size' : effpot_fontsize,
                                            }
                                        ),                                        
                                        html.Br(),
                                        html.Label("Periapsis: ",
                                                   style={
                                                       'font-size': effpot_fontsize,
                                                       'margin-top': '10px',
                                                       'margin-right': 10,
                                                       'margin-bottom': '10px',
                                                   }),
                                        dcc.Input(
                                            id='periap-val-str', type='text',
                                            value = default_periap_str,
                                            placeholder = default_periap_str,
                                            # don't allow update until user enters
                                            debounce = True,
                                            size = '7',
                                            style={
                                                'display': 'inline-block',
                                                'margin-right': 2,
                                                'font-size' : effpot_fontsize,
                                            }
                                        ),
                                        html.Div(
                                            children=" > 4 M",
                                            style={
                                                'display': 'inline-block',
                                                'margin-left': 5,
                                                'font-size' : effpot_fontsize,
                                            }
                                        ),
                                    ]),
                                    html.P("Units are \"geometrized\" (G=c=1)."
                                           + " Only bound orbits are allowed."
                                           + " Upper limit on angular momentum"
                                           + " is for computational efficiency; all"
                                           + " other limits are physical.",
                                        style={
                                            #'display': 'inline-block',
                                            'font-size' : effpot_fontsize_small,
                                            'margin-top' : '10px',
                                            'margin-bottom' :'5px'
                                        }
                                    ),
                                ]),
                            ]),
                            dbc.Col(width=6, children=[
                                dcc.Graph(
                                    id='orbitgraph',
                                    figure=init_orbit_fig,
                                    style={
                                        'height': '95%',
                                        'margin-right': '10px',
                                        'margin-left': '10px',
                                        #'float': 'right'
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                )
                            ]),
                        ]),
                    ]),
                    dbc.Card(
                        [
                            dbc.CardHeader("Effective Potential"),
                            dbc.CardBody([
                                dcc.Graph(
                                    id='potentialgraph',
                                    figure=init_pot_fig,
                                    style={
                                        'height': '75%',
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                                html.Div( style={
                                    'textAlign': 'center',
                                    'margin-top': '10px',
                                    'border': '2px lightgrey solid',
                                                 }, children=[
                                    html.Button('revert to default', id='default-button', n_clicks=0,
                                                style ={
                                                    'font-size' : effpot_fontsize_small,
                                                }),
                                    html.Div(
                                        children="Angular Momentum:",
                                        style={
                                            'display': 'inline-block',
                                            'margin-left': 15,
                                            'margin-right': 10,
                                            'margin-top': 10,
                                            'margin-bottom': 10,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    html.Div(
                                        children= "3.464 M < ",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 5,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    dcc.Input(
                                        id='angmom-val-str', type='text',
                                        value = default_angmom_str,
                                        placeholder = default_angmom_str,
                                        # don't allow update until user enters
                                        debounce = True,
                                        size = '7',
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 5,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    html.Div(
                                        children=" < " + ell_max_str + " M",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 30,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    html.Div(
                                        children="Energy: ",
                                        style={
                                            'display': 'inline-block',
                                            'margin-right': 10,
                                            'font-size' : effpot_fontsize,
                                        }
                                    ),
                                    dcc.Input(
                                        id='energy-val-str', type='text',
                                        value = default_energy_str,
                                        placeholder = default_energy_str,
                                        # don't allow update until user enters
                                        debounce = True,
                                        size = '7',
                                        style={
                                            'display': 'inline-block',
                                            'font-size' : effpot_fontsize,
                                            'margin-right': 5,
                                        }
                                    ),
                                ]),
                                html.Div(children="Hint: First try adjusting energy; adjusting angular momentum changes the Veff function (and the view).",
                                        style={
                                            'font-size' : effpot_fontsize_small,
                                            'font-weight': 'bold',
                                            'color': '#636EFA',
                                            'margin-top' : '10px',
                                        }
                                    ),                              
                                
                            ]),
                        ],
                    )
                    ],
                )
            ]),
            #
            #=== Second Row
            #
            dbc.Row([
                dbc.Card(
                    [
                        dbc.CardHeader([
                            html.Div(children="Gravitational Wave Signal for M/m = ",
                                     style={'display': 'inline-block'}),
                            dcc.Input(
                                id='bhmass-str', type='text',
                                value = default_bhmass_str,
                                placeholder = default_bhmass_str,
                                # don't allow update until user enters
                                debounce = True,
                                size = '8',
                                style={
                                    'display': 'inline-block',
                                    'margin-right': 2,
                                    'font-size' : effpot_fontsize,
                                }
                            ),
                            html.Div(children=" \u2265 " + bhmass_min_str,  # \u2265 is ≥
                                     style={'display': 'inline-block'}),
                        ]),
                        dbc.CardBody(
                            [
                                dcc.Graph(
                                    id='gw_plus_graph',
                                    figure=init_gw_plus_fig,
                                    style={
                                        'display': 'inline-block',
                                        'width' : '50%',
                                        # should be able to say height=100%
                                        # but it does not work, so this must be
                                        # something smaller than the card height
                                        # specified below in "vh" (viewwindow height)
                                        'height' : height_secondrow,
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                                dcc.Graph(
                                    id='gw_cross_graph',
                                    figure=init_gw_cross_fig,
                                    style={
                                        'display': 'inline-block',
                                        'width' : '50%',
                                        'height' : height_secondrow,
                                    },
                                    config={
                                        'displayModeBar': False,
                                    }
                                ),
                            ]
                        )
                    ],
                    style = { 'height' : height_secondrow_card }
                )
            ])
        ],
    )
    return dashboard_page
###############################
# end dashboard webpage setup #
###############################
