/************************************************************
 *** Functions creating scripts for plotting with gnuplot ***
 ************************************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

 /*** writePalViridis ***
  * Writes a palette file for gnuplot ("viridis" from Matplotlib) */
void writePalViridis()
{
	ofstream pal("palette.pal");
	pal << "# Viridis color map" << endl << endl;
	pal << "set palette defined (\\" << endl;
	pal << "0	0.267004   0.004874   0.329415,\\" << endl;
	pal << "1	0.268510   0.009605   0.335427,\\" << endl;
	pal << "2	0.269944   0.014625   0.341379,\\" << endl;
	pal << "3	0.271305   0.019942   0.347269,\\" << endl;
	pal << "4	0.272594   0.025563   0.353093,\\" << endl;
	pal << "5	0.273809   0.031497   0.358853,\\" << endl;
	pal << "6	0.274952   0.037752   0.364543,\\" << endl;
	pal << "7	0.276022   0.044167   0.370164,\\" << endl;
	pal << "8	0.277018   0.050344   0.375715,\\" << endl;
	pal << "9	0.277941   0.056324   0.381191,\\" << endl;
	pal << "10	0.278791   0.062145   0.386592,\\" << endl;
	pal << "11	0.279566   0.067836   0.391917,\\" << endl;
	pal << "12	0.280267   0.073417   0.397163,\\" << endl;
	pal << "13	0.280894   0.078907   0.402329,\\" << endl;
	pal << "14	0.281446   0.084320   0.407414,\\" << endl;
	pal << "15	0.281924   0.089666   0.412415,\\" << endl;
	pal << "16	0.282327   0.094955   0.417331,\\" << endl;
	pal << "17	0.282656   0.100196   0.422160,\\" << endl;
	pal << "18	0.282910   0.105393   0.426902,\\" << endl;
	pal << "19	0.283091   0.110553   0.431554,\\" << endl;
	pal << "20	0.283197   0.115680   0.436115,\\" << endl;
	pal << "21	0.283229   0.120777   0.440584,\\" << endl;
	pal << "22	0.283187   0.125848   0.444960,\\" << endl;
	pal << "23	0.283072   0.130895   0.449241,\\" << endl;
	pal << "24	0.282884   0.135920   0.453427,\\" << endl;
	pal << "25	0.282623   0.140926   0.457517,\\" << endl;
	pal << "26	0.282290   0.145912   0.461510,\\" << endl;
	pal << "27	0.281887   0.150881   0.465405,\\" << endl;
	pal << "28	0.281412   0.155834   0.469201,\\" << endl;
	pal << "29	0.280868   0.160771   0.472899,\\" << endl;
	pal << "30	0.280255   0.165693   0.476498,\\" << endl;
	pal << "31	0.279574   0.170599   0.479997,\\" << endl;
	pal << "32	0.278826   0.175490   0.483397,\\" << endl;
	pal << "33	0.278012   0.180367   0.486697,\\" << endl;
	pal << "34	0.277134   0.185228   0.489898,\\" << endl;
	pal << "35	0.276194   0.190074   0.493001,\\" << endl;
	pal << "36	0.275191   0.194905   0.496005,\\" << endl;
	pal << "37	0.274128   0.199721   0.498911,\\" << endl;
	pal << "38	0.273006   0.204520   0.501721,\\" << endl;
	pal << "39	0.271828   0.209303   0.504434,\\" << endl;
	pal << "40	0.270595   0.214069   0.507052,\\" << endl;
	pal << "41	0.269308   0.218818   0.509577,\\" << endl;
	pal << "42	0.267968   0.223549   0.512008,\\" << endl;
	pal << "43	0.266580   0.228262   0.514349,\\" << endl;
	pal << "44	0.265145   0.232956   0.516599,\\" << endl;
	pal << "45	0.263663   0.237631   0.518762,\\" << endl;
	pal << "46	0.262138   0.242286   0.520837,\\" << endl;
	pal << "47	0.260571   0.246922   0.522828,\\" << endl;
	pal << "48	0.258965   0.251537   0.524736,\\" << endl;
	pal << "49	0.257322   0.256130   0.526563,\\" << endl;
	pal << "50	0.255645   0.260703   0.528312,\\" << endl;
	pal << "51	0.253935   0.265254   0.529983,\\" << endl;
	pal << "52	0.252194   0.269783   0.531579,\\" << endl;
	pal << "53	0.250425   0.274290   0.533103,\\" << endl;
	pal << "54	0.248629   0.278775   0.534556,\\" << endl;
	pal << "55	0.246811   0.283237   0.535941,\\" << endl;
	pal << "56	0.244972   0.287675   0.537260,\\" << endl;
	pal << "57	0.243113   0.292092   0.538516,\\" << endl;
	pal << "58	0.241237   0.296485   0.539709,\\" << endl;
	pal << "59	0.239346   0.300855   0.540844,\\" << endl;
	pal << "60	0.237441   0.305202   0.541921,\\" << endl;
	pal << "61	0.235526   0.309527   0.542944,\\" << endl;
	pal << "62	0.233603   0.313828   0.543914,\\" << endl;
	pal << "63	0.231674   0.318106   0.544834,\\" << endl;
	pal << "64	0.229739   0.322361   0.545706,\\" << endl;
	pal << "65	0.227802   0.326594   0.546532,\\" << endl;
	pal << "66	0.225863   0.330805   0.547314,\\" << endl;
	pal << "67	0.223925   0.334994   0.548053,\\" << endl;
	pal << "68	0.221989   0.339161   0.548752,\\" << endl;
	pal << "69	0.220057   0.343307   0.549413,\\" << endl;
	pal << "70	0.218130   0.347432   0.550038,\\" << endl;
	pal << "71	0.216210   0.351535   0.550627,\\" << endl;
	pal << "72	0.214298   0.355619   0.551184,\\" << endl;
	pal << "73	0.212395   0.359683   0.551710,\\" << endl;
	pal << "74	0.210503   0.363727   0.552206,\\" << endl;
	pal << "75	0.208623   0.367752   0.552675,\\" << endl;
	pal << "76	0.206756   0.371758   0.553117,\\" << endl;
	pal << "77	0.204903   0.375746   0.553533,\\" << endl;
	pal << "78	0.203063   0.379716   0.553925,\\" << endl;
	pal << "79	0.201239   0.383670   0.554294,\\" << endl;
	pal << "80	0.199430   0.387607   0.554642,\\" << endl;
	pal << "81	0.197636   0.391528   0.554969,\\" << endl;
	pal << "82	0.195860   0.395433   0.555276,\\" << endl;
	pal << "83	0.194100   0.399323   0.555565,\\" << endl;
	pal << "84	0.192357   0.403199   0.555836,\\" << endl;
	pal << "85	0.190631   0.407061   0.556089,\\" << endl;
	pal << "86	0.188923   0.410910   0.556326,\\" << endl;
	pal << "87	0.187231   0.414746   0.556547,\\" << endl;
	pal << "88	0.185556   0.418570   0.556753,\\" << endl;
	pal << "89	0.183898   0.422383   0.556944,\\" << endl;
	pal << "90	0.182256   0.426184   0.557120,\\" << endl;
	pal << "91	0.180629   0.429975   0.557282,\\" << endl;
	pal << "92	0.179019   0.433756   0.557430,\\" << endl;
	pal << "93	0.177423   0.437527   0.557565,\\" << endl;
	pal << "94	0.175841   0.441290   0.557685,\\" << endl;
	pal << "95	0.174274   0.445044   0.557792,\\" << endl;
	pal << "96	0.172719   0.448791   0.557885,\\" << endl;
	pal << "97	0.171176   0.452530   0.557965,\\" << endl;
	pal << "98	0.169646   0.456262   0.558030,\\" << endl;
	pal << "99	0.168126   0.459988   0.558082,\\" << endl;
	pal << "100	0.166617   0.463708   0.558119,\\" << endl;
	pal << "101	0.165117   0.467423   0.558141,\\" << endl;
	pal << "102	0.163625   0.471133   0.558148,\\" << endl;
	pal << "103	0.162142   0.474838   0.558140,\\" << endl;
	pal << "104	0.160665   0.478540   0.558115,\\" << endl;
	pal << "105	0.159194   0.482237   0.558073,\\" << endl;
	pal << "106	0.157729   0.485932   0.558013,\\" << endl;
	pal << "107	0.156270   0.489624   0.557936,\\" << endl;
	pal << "108	0.154815   0.493313   0.557840,\\" << endl;
	pal << "109	0.153364   0.497000   0.557724,\\" << endl;
	pal << "110	0.151918   0.500685   0.557587,\\" << endl;
	pal << "111	0.150476   0.504369   0.557430,\\" << endl;
	pal << "112	0.149039   0.508051   0.557250,\\" << endl;
	pal << "113	0.147607   0.511733   0.557049,\\" << endl;
	pal << "114	0.146180   0.515413   0.556823,\\" << endl;
	pal << "115	0.144759   0.519093   0.556572,\\" << endl;
	pal << "116	0.143343   0.522773   0.556295,\\" << endl;
	pal << "117	0.141935   0.526453   0.555991,\\" << endl;
	pal << "118	0.140536   0.530132   0.555659,\\" << endl;
	pal << "119	0.139147   0.533812   0.555298,\\" << endl;
	pal << "120	0.137770   0.537492   0.554906,\\" << endl;
	pal << "121	0.136408   0.541173   0.554483,\\" << endl;
	pal << "122	0.135066   0.544853   0.554029,\\" << endl;
	pal << "123	0.133743   0.548535   0.553541,\\" << endl;
	pal << "124	0.132444   0.552216   0.553018,\\" << endl;
	pal << "125	0.131172   0.555899   0.552459,\\" << endl;
	pal << "126	0.129933   0.559582   0.551864,\\" << endl;
	pal << "127	0.128729   0.563265   0.551229,\\" << endl;
	pal << "128	0.127568   0.566949   0.550556,\\" << endl;
	pal << "129	0.126453   0.570633   0.549841,\\" << endl;
	pal << "130	0.125394   0.574318   0.549086,\\" << endl;
	pal << "131	0.124395   0.578002   0.548287,\\" << endl;
	pal << "132	0.123463   0.581687   0.547445,\\" << endl;
	pal << "133	0.122606   0.585371   0.546557,\\" << endl;
	pal << "134	0.121831   0.589055   0.545623,\\" << endl;
	pal << "135	0.121148   0.592739   0.544641,\\" << endl;
	pal << "136	0.120565   0.596422   0.543611,\\" << endl;
	pal << "137	0.120092   0.600104   0.542530,\\" << endl;
	pal << "138	0.119738   0.603785   0.541400,\\" << endl;
	pal << "139	0.119512   0.607464   0.540218,\\" << endl;
	pal << "140	0.119423   0.611141   0.538982,\\" << endl;
	pal << "141	0.119483   0.614817   0.537692,\\" << endl;
	pal << "142	0.119699   0.618490   0.536347,\\" << endl;
	pal << "143	0.120081   0.622161   0.534946,\\" << endl;
	pal << "144	0.120638   0.625828   0.533488,\\" << endl;
	pal << "145	0.121380   0.629492   0.531973,\\" << endl;
	pal << "146	0.122312   0.633153   0.530398,\\" << endl;
	pal << "147	0.123444   0.636809   0.528763,\\" << endl;
	pal << "148	0.124780   0.640461   0.527068,\\" << endl;
	pal << "149	0.126326   0.644107   0.525311,\\" << endl;
	pal << "150	0.128087   0.647749   0.523491,\\" << endl;
	pal << "151	0.130067   0.651384   0.521608,\\" << endl;
	pal << "152	0.132268   0.655014   0.519661,\\" << endl;
	pal << "153	0.134692   0.658636   0.517649,\\" << endl;
	pal << "154	0.137339   0.662252   0.515571,\\" << endl;
	pal << "155	0.140210   0.665859   0.513427,\\" << endl;
	pal << "156	0.143303   0.669459   0.511215,\\" << endl;
	pal << "157	0.146616   0.673050   0.508936,\\" << endl;
	pal << "158	0.150148   0.676631   0.506589,\\" << endl;
	pal << "159	0.153894   0.680203   0.504172,\\" << endl;
	pal << "160	0.157851   0.683765   0.501686,\\" << endl;
	pal << "161	0.162016   0.687316   0.499129,\\" << endl;
	pal << "162	0.166383   0.690856   0.496502,\\" << endl;
	pal << "163	0.170948   0.694384   0.493803,\\" << endl;
	pal << "164	0.175707   0.697900   0.491033,\\" << endl;
	pal << "165	0.180653   0.701402   0.488189,\\" << endl;
	pal << "166	0.185783   0.704891   0.485273,\\" << endl;
	pal << "167	0.191090   0.708366   0.482284,\\" << endl;
	pal << "168	0.196571   0.711827   0.479221,\\" << endl;
	pal << "169	0.202219   0.715272   0.476084,\\" << endl;
	pal << "170	0.208030   0.718701   0.472873,\\" << endl;
	pal << "171	0.214000   0.722114   0.469588,\\" << endl;
	pal << "172	0.220124   0.725509   0.466226,\\" << endl;
	pal << "173	0.226397   0.728888   0.462789,\\" << endl;
	pal << "174	0.232815   0.732247   0.459277,\\" << endl;
	pal << "175	0.239374   0.735588   0.455688,\\" << endl;
	pal << "176	0.246070   0.738910   0.452024,\\" << endl;
	pal << "177	0.252899   0.742211   0.448284,\\" << endl;
	pal << "178	0.259857   0.745492   0.444467,\\" << endl;
	pal << "179	0.266941   0.748751   0.440573,\\" << endl;
	pal << "180	0.274149   0.751988   0.436601,\\" << endl;
	pal << "181	0.281477   0.755203   0.432552,\\" << endl;
	pal << "182	0.288921   0.758394   0.428426,\\" << endl;
	pal << "183	0.296479   0.761561   0.424223,\\" << endl;
	pal << "184	0.304148   0.764704   0.419943,\\" << endl;
	pal << "185	0.311925   0.767822   0.415586,\\" << endl;
	pal << "186	0.319809   0.770914   0.411152,\\" << endl;
	pal << "187	0.327796   0.773980   0.406640,\\" << endl;
	pal << "188	0.335885   0.777018   0.402049,\\" << endl;
	pal << "189	0.344074   0.780029   0.397381,\\" << endl;
	pal << "190	0.352360   0.783011   0.392636,\\" << endl;
	pal << "191	0.360741   0.785964   0.387814,\\" << endl;
	pal << "192	0.369214   0.788888   0.382914,\\" << endl;
	pal << "193	0.377779   0.791781   0.377939,\\" << endl;
	pal << "194	0.386433   0.794644   0.372886,\\" << endl;
	pal << "195	0.395174   0.797475   0.367757,\\" << endl;
	pal << "196	0.404001   0.800275   0.362552,\\" << endl;
	pal << "197	0.412913   0.803041   0.357269,\\" << endl;
	pal << "198	0.421908   0.805774   0.351910,\\" << endl;
	pal << "199	0.430983   0.808473   0.346476,\\" << endl;
	pal << "200	0.440137   0.811138   0.340967,\\" << endl;
	pal << "201	0.449368   0.813768   0.335384,\\" << endl;
	pal << "202	0.458674   0.816363   0.329727,\\" << endl;
	pal << "203	0.468053   0.818921   0.323998,\\" << endl;
	pal << "204	0.477504   0.821444   0.318195,\\" << endl;
	pal << "205	0.487026   0.823929   0.312321,\\" << endl;
	pal << "206	0.496615   0.826376   0.306377,\\" << endl;
	pal << "207	0.506271   0.828786   0.300362,\\" << endl;
	pal << "208	0.515992   0.831158   0.294279,\\" << endl;
	pal << "209	0.525776   0.833491   0.288127,\\" << endl;
	pal << "210	0.535621   0.835785   0.281908,\\" << endl;
	pal << "211	0.545524   0.838039   0.275626,\\" << endl;
	pal << "212	0.555484   0.840254   0.269281,\\" << endl;
	pal << "213	0.565498   0.842430   0.262877,\\" << endl;
	pal << "214	0.575563   0.844566   0.256415,\\" << endl;
	pal << "215	0.585678   0.846661   0.249897,\\" << endl;
	pal << "216	0.595839   0.848717   0.243329,\\" << endl;
	pal << "217	0.606045   0.850733   0.236712,\\" << endl;
	pal << "218	0.616293   0.852709   0.230052,\\" << endl;
	pal << "219	0.626579   0.854645   0.223353,\\" << endl;
	pal << "220	0.636902   0.856542   0.216620,\\" << endl;
	pal << "221	0.647257   0.858400   0.209861,\\" << endl;
	pal << "222	0.657642   0.860219   0.203082,\\" << endl;
	pal << "223	0.668054   0.861999   0.196293,\\" << endl;
	pal << "224	0.678489   0.863742   0.189503,\\" << endl;
	pal << "225	0.688944   0.865448   0.182725,\\" << endl;
	pal << "226	0.699415   0.867117   0.175971,\\" << endl;
	pal << "227	0.709898   0.868751   0.169257,\\" << endl;
	pal << "228	0.720391   0.870350   0.162603,\\" << endl;
	pal << "229	0.730889   0.871916   0.156029,\\" << endl;
	pal << "230	0.741388   0.873449   0.149561,\\" << endl;
	pal << "231	0.751884   0.874951   0.143228,\\" << endl;
	pal << "232	0.762373   0.876424   0.137064,\\" << endl;
	pal << "233	0.772852   0.877868   0.131109,\\" << endl;
	pal << "234	0.783315   0.879285   0.125405,\\" << endl;
	pal << "235	0.793760   0.880678   0.120005,\\" << endl;
	pal << "236	0.804182   0.882046   0.114965,\\" << endl;
	pal << "237	0.814576   0.883393   0.110347,\\" << endl;
	pal << "238	0.824940   0.884720   0.106217,\\" << endl;
	pal << "239	0.835270   0.886029   0.102646,\\" << endl;
	pal << "240	0.845561   0.887322   0.099702,\\" << endl;
	pal << "241	0.855810   0.888601   0.097452,\\" << endl;
	pal << "242	0.866013   0.889868   0.095953,\\" << endl;
	pal << "243	0.876168   0.891125   0.095250,\\" << endl;
	pal << "244	0.886271   0.892374   0.095374,\\" << endl;
	pal << "245	0.896320   0.893616   0.096335,\\" << endl;
	pal << "246	0.906311   0.894855   0.098125,\\" << endl;
	pal << "247	0.916242   0.896091   0.100717,\\" << endl;
	pal << "248	0.926106   0.897330   0.104071,\\" << endl;
	pal << "249	0.935904   0.898570   0.108131,\\" << endl;
	pal << "250	0.945636   0.899815   0.112838,\\" << endl;
	pal << "251	0.955300   0.901065   0.118128,\\" << endl;
	pal << "252	0.964894   0.902323   0.123941,\\" << endl;
	pal << "253	0.974417   0.903590   0.130215,\\" << endl;
	pal << "254	0.983868   0.904867   0.136897,\\" << endl;
	pal << "255	0.993248   0.906157   0.143936)" << endl;

	pal.close();
}


/*** createNetworkColorPlot ***
 * Creates a network color plot of data that has to be specified *
 * - f: the gnuplot file to write to
 * - Nl: the number of neurons in one line
 * - irradiance: the irradiance (amplitude) used for the simulation (set this to a negative value if irrelevant)
 * - column: the column in the data file to be used for color-coded values
 * - prefix: the principal name of the plot
 * - postfix: if necessary, the type of stimulation (actually either "_light" or "_current")
 * - matrix: if true, a non-interpolated matrix plot is created, else an interpolated color plot
 * - cblabel: the label for the colorbar
 * - cbmin [optional]: the minimum for the colorbar
 * - cbmax [optional]: the maximum for the colorbar
 * - cbtics [optional]: one increment for the colorbar axis */
void createNetworkColorPlot(ofstream &f, int Nl, double irradiance, int column, const char* prefix, const char* postfix, bool matrix, const char* cblabel,
							double cbmin = -1.0, double cbmax = -1.0, double cbtics = -1.0)
{
	char gplcall[100];

	if (irradiance >= 0.0)
		sprintf(gplcall, "%s_%.2f", prefix, irradiance);
	else
		sprintf(gplcall, "%s", prefix);

	f << "# Output configuration" << endl;
	f << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5 size 5,4.5" << endl;
	f << "set output '" << dateStr("_") << gplcall << "_map" << postfix << ".pdf'" << endl << endl;

	f << "# Contour plot configuration" << endl;
	if (!matrix)
	{
		f << "set pm3d" << endl;
		f << "unset surface" << endl;
	}
	f << "set view map	# 2-dim. plot" << endl;
	f << "#set view 50,75 # just for 3-dim. plot to rotate # 3-dim. plot" << endl;
	f << "#set ticslevel 0 # just for 3-dim. plot to shift z=0 down into the xy plane # 3-dim. plot" << endl;
	f << "unset key" << endl;
	if (!matrix)
		f << "set pm3d interpolate 100,100 # interpolate the color (0: gnuplot automatically decides the steps to use)" << endl;
	f << "load 'palette.pal' # load palette" << endl;
	f.precision(5); // set number of significant digits in output
	if (strcmp(prefix, "firingrate") == 0)
		f << "set title \"Firing rates across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "irradiance") == 0)
		f << "set title \"Light intensities across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "inc_connection") == 0)
		f << "set title \"Incoming connections per neuron" << endl;
	f << "set size square" << endl << endl;  // quadratic shape

	f << "# Axis configuration" << endl;
	if (!matrix)
	{
		f << "set xrange [1" << ":" << Nl << "]" << endl;
		f << "set yrange [1" << ":" << Nl << "]" << endl;
	}
	else
	{
		f << "set xrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
		f << "set yrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
	}
	f << "set mxtics 2" << endl;
	f << "set mytics 2" << endl;
	if (cbmin >= 0.0 && cbmax >= 0.0)
	{
		if (cbmin == cbmax)
			f << "set cbrange [" << 0.95*cbmin << ":" << 1.05*cbmin << "]" << endl;
		else
			f << "set cbrange [" << cbmin << ":" << cbmax << "]" << endl;
	}
	if (cbtics > 0.0)
		f << "set cbtics " << cbtics << endl;
	f << "set xlabel \"Neuron column\"" << endl;
	f << "set ylabel \"Neuron row\" offset -1.5" << endl; // offset to shift the label to the left
	f << "set cblabel \"" << cblabel << "\"" << endl;
	f << "set format x '%1.0f'" << endl;
	f << "set format y '%1.0f'" << endl;
	f << "set format z '%.5f'" << endl << endl;

	f << "# Set margins" << endl;
	f << "set tmargin at screen 0.90" << endl;
	f << "set bmargin at screen 0.17" << endl;
	f << "set rmargin at screen 0.78" << endl;
	f << "set lmargin at screen 0.15" << endl << endl;

	f << "# Plotting" << endl;
	f << "set datafile missing \"nan\"" << endl;
	if (!matrix)
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " notitle with lines lt 1" << endl;
	else
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " with image" << endl;
	f.close();

	if (irradiance >= 0.0)
		sprintf(gplcall, "gnuplot %s_%.2f_map%s.gpl", prefix, irradiance, postfix);
	else
		sprintf(gplcall, "gnuplot %s_map%s.gpl", prefix, postfix);
	system(gplcall);
}


/*** createSpikeNumberPlot ***
 * Creates a plot of the spike number in a time bin over time *
 * - gpl: handle to a file containing all the gnuplot calls *
 * - t_max: the duration of the simulation in s *
 * - bin: the size of one time bin in s */
void createSpikeNumberPlot(ofstream &gpl, double t_max, double bin)
{
	ofstream f("spike_number.gpl");
	if (!f.is_open())
		return;

	f << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5" << endl
	  << "set output '" << dateStr("_spike_number.pdf'") << endl
	  << "set xlabel \"time / s\"" << endl
	  << "set ylabel \"spikes per second\"" << endl
	  << "set tmargin at screen 0.95" << endl
	  << "set bmargin at screen 0.23" << endl << endl
	  << "plot [x=0:" << t_max << "] '" << dateStr("_spike_number.txt") << "' \\" << endl
	  << "     using 1:($2/" << bin << ") notitle with lines";

	// Call gnuplot script
	f.close();
	string gplc = "gnuplot spike_number.gpl";
	gpl << gplc << endl;
	system(gplc.c_str());
}

/*** createSpikeRasterPlot ***
 * Creates a spike raster plot of given data and offers the possibility to emphasize some neurons by green color *
 * - gpl: handle to a file containing all the gnuplot calls *
 * - t_start: the time minimum for plotting *
 * - t_end: the time maximum for plotting *
 * - Ne: the number of excitatory neurons *
 * - Ni: the number of inhibitory neurons *
 * - Nc_start: the number of the first emphasized neuron *
 * - Nc_end: the number of the last emphasized neuron plus one */
void createSpikeRasterPlot(ofstream &gpl, double t_start, double t_end, double Ne, double Ni, double Nc_start, double Nc_end)
{
	int num_ytics; // number of tics at y-axis
	int yn_step; // distance in neurons between tics at y-axis
	double ypos = 0.05; // start, bottom and top margin size
	double ypos_step; // go in steps up to y=1.0-ypos

	ofstream f("spike_raster.gpl");

	if (!f.is_open())
		return;

	if (Ne+Ni > 5)
		num_ytics = 5;
	else
		num_ytics = Ne+Ni;

	yn_step = int((Ne+Ni)/double(num_ytics));
	ypos_step = (1.0-2*ypos)/num_ytics;

	// Set PDF output
	f << "#set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;
	f << "#set output '" << dateStr("_spike_raster.pdf'") << endl << endl;

	// Set PNG output
	f << "set term png enhanced font Sans 20 size 1280,960 lw 2.5" << endl;
	f << "set output '" << dateStr("_spike_raster.png'") << endl << endl;

	// Assign variables for neuron number
	f << "Ne = " << Ne << endl;
	f << "Ni = " << Ni << endl;
	f << "Nc_start = " << Nc_start << endl;
	f << "Nc_end = " << Nc_end << endl;

	// Set labels
	f << "set xlabel \"t / s\"" << endl;
	f << "unset ylabel" << endl;
	f << "set yrange [0:1]" << endl;
	f << "set ytics out (";

	for (int i=0; i<Ne+Ni; i++)
	{
		if (i % yn_step == 0)
		{
			int tic_number = i / yn_step;
			f << "\"#" << i << "\" " << dtos(ypos+tic_number*ypos_step, 2);
			if (tic_number + 1 < num_ytics) // if the last y-tic is not yet reached
				f << ", ";
		}
	}
	f << ")" << endl << endl;

	// Plot spike time series
	f << "plot [x=" << t_start << ":" << t_end << "] '" << dateStr("_spike_raster.txt") << "' using 1:($2 < Ne ? (" << 1.-2*ypos << "*$2/(Ne+Ni) + " << ypos << ") : 1/0) "
	  //<< "notitle with points pt 7 ps 0.01 lc \"blue\", \\" << endl
	  << "notitle with dots lc \"blue\", \\" << endl
	  << "     '" << dateStr("_spike_raster.txt") << "' using 1:($2 >= Ne ? (" << 1.-2*ypos << "*$2/(Ne+Ni) + " << ypos << ") : 1/0) "
	  << "notitle with dots lc \"red\"";
	  //<< ", \\" << endl << "     '" << dateStr("_spike_raster.txt") << "' using 1:(($2 >= Nc_start && $2 < Nc_end) ? (" << 1.-2*ypos << "*$2/(Ne+Ni) + " << ypos << ") : 1/0) "
	  //<< "notitle with dots lc \"green\""; // additional emphasis on the na subpopulation

	// Call gnuplot script
	f.close();
	string gplc = "gnuplot spike_raster.gpl";
	gpl << gplc << endl;
	system(gplc.c_str());
}


/*** createMeanWeightPlotCA ***
 * Creates a plot that displays the mean early- and late-phase and the total strength as well as the protein amount in the cell assembly over time *
 * - gpl: handle to a file containing the gnuplot calls *
 * - t_max: the duration of the time series in s *
 * - h_0: initial weight */
void createMeanWeightPlotCA(ofstream &gpl, double t_max, double h_0)
{
	int x_max; // maximum x-value
	double x_scale = 100./96.; // value to enlarge x-axis
	string x_data; // string containing the column specifier for x-values
	string x_label = "Time "; // string containing the x-label

	if (t_max > 3600)
	{
		x_max = round(x_scale*t_max/60.);
		x_data = "($1/60)";
		x_label += "(min)";
	}
	else
	{
		x_max = round(x_scale*t_max);
		x_data = "1";
		x_label += "(s)";
	}

	string spaces = "                                            ";
	string filename = string("mean_weight_plot_CA.gpl");
	ofstream f(filename);
	if (!f.is_open())
		return;

	f << "set term pdf enhanced font \"Sans, 18\"" << endl
	  << "set output \"" << dateStr("_mean_weight_CA.pdf\"") << endl
	  << endl
	  << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl
	  << endl
	  << "unset key" << endl // for better view
	  << "set ytics 1 nomirror" << endl
	  << endl
	  << "h0 = " << h_0 << endl
	  << endl
	  << "set xlabel \"" << x_label << "\"" << endl
	  << "set ylabel \"arb. units\"" << endl
	  << endl
		<< "set style fill transparent solid 0.8 noborder" << endl
	  << "set samples 1000" << endl
	  << "plot [0:" << x_max << "][0.0:3.2] "
		<< "\"" << dateStr("_mean_weight.txt") << "\" using " << x_data << ":(($2+h0*$4-$3-h0*$5)/h0):(($2+h0*$4+$3+h0*$5)/h0) notitle lw 2 lc rgb rgb(255,127,14) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":(($2-$3)/h0):(($2+$3)/h0) notitle lw 2 lc rgb rgb(128,0,0) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($4-$5+1):($4+$5+1) notitle lw 2 lc rgb rgb(31,119,180) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($6-$7):($6+$7) notitle lw 2 lc rgb rgb (0,128,0) with filledcurves, \\" << endl << spaces
                << "\"\" using " << x_data << ":(($2+h0*$4)/h0) title \"⟨w(t)/h_0⟩\" lw 2 lc rgb rgb(255,127,14) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($2/h0) title \"⟨h(t)/h_0⟩\" lw 2 lc rgb rgb(128,0,0) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($4+1) title \"⟨z(t)+1⟩\" lw 2 lc rgb rgb(31,119,180) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":6 title \"⟨p(t)⟩\" lw 2 lc rgb rgb(0,128,0) with lines, \\" << endl << spaces;


	// Call gnuplot script
	f.close();
	string gplc = concat("gnuplot ", filename);
	gpl << gplc << endl;
	system(gplc.c_str());
}

/*** createMeanWeightPlotControl ***
 * Creates a plot that displays the mean early- and late-phase and the total strength as well as the protein amount in the control subpopulation over time *
 * - gpl: handle to a file containing the gnuplot calls *
 * - t_max: the duration of the time series in s *
 * - h_0: initial weight */
void createMeanWeightPlotControl(ofstream &gpl, double t_max, double h_0)
{
	int x_max; // maximum x-value
	double x_scale = 100./96.; // value to enlarge x-axis
	string x_data; // string containing the column specifier for x-values
	string x_label = "Time "; // string containing the x-label

	if (t_max > 3600)
	{
		x_max = round(x_scale*t_max/60.);
		x_data = "($1/60)";
		x_label += "(min)";
	}
	else
	{
		x_max = round(x_scale*t_max);
		x_data = "1";
		x_label += "(s)";
	}

	string spaces = "                                            ";
	string filename = string("mean_weight_plot_control.gpl");
	ofstream f(filename);
	if (!f.is_open())
		return;

	f << "set term pdf enhanced font \"Sans, 18\"" << endl
	  << "set output \"" << dateStr("_mean_weight_control.pdf\"") << endl
	  << endl
	  << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl
	  << endl
	  << "unset key" << endl // for better view
	  << "set ytics 1 nomirror" << endl
	  << endl
	  << "h0 = " << h_0 << endl
	  << endl
	  << "set xlabel \"" << x_label << "\"" << endl
	  << "set ylabel \"arb. units\"" << endl
	  << endl
		<< "set style fill transparent solid 0.8 noborder" << endl
	  << "set samples 1000" << endl
	  << "plot [0:" << x_max << "][0.0:3.2] "
		<< "\"" << dateStr("_mean_weight.txt") << "\" using " << x_data << ":(($8+h0*$10-$9-h0*$11)/h0):(($8+h0*$10+$9+h0*$11)/h0) notitle lw 2 lc rgb rgb(255,127,14) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":(($8-$9)/h0):(($8+$9)/h0) notitle lw 2 lc rgb rgb(128,0,0) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($10-$11+1):($10+$11+1) notitle lw 2 lc rgb rgb(31,119,180) with filledcurves, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($12-$13):($12+$13) notitle lw 2 lc rgb rgb (0,128,0) with filledcurves, \\" << endl << spaces
                << "\"\" using " << x_data << ":(($8+h0*$10)/h0) title \"⟨w(t)/h_0⟩\" lw 2 lc rgb rgb(255,127,14) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($8/h0) title \"⟨h(t)/h_0⟩\" lw 2 lc rgb rgb(128,0,0) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":($10+1) title \"⟨z(t)+1⟩\" lw 2 lc rgb rgb(31,119,180) with lines, \\" << endl << spaces
	        << "\"\" using " << x_data << ":12 title \"⟨p(t)⟩\" lw 2 lc rgb rgb(0,128,0) with lines, \\" << endl << spaces;


	// Call gnuplot script
	f.close();
	string gplc = concat("gnuplot ", filename);
	gpl << gplc << endl;
	system(gplc.c_str());
}


/*** createSynapsePlot ***
 * Creates a plot that displays synapse-specific data (early-/late-phase plasticity, total strength) over time *
 * - synapse_output: a vector of the considered synapses *
 * - gpl: handle to a file containing the gnuplot calls *
 * - t_max: the duration of the time series in s *
 * - stimulus: a label for the used stimulus protocol *
 * - h_0: initial weight *
 * - theta_tag: tagging threshold *
 * - theta_pro: protein synthesis threshold *
 * - theta_p: potentiation threshold for Ca dynamics *
 * - theta_d: depression threshold for Ca dynamics *
 * - plh_cols: number of columns to skip for they contain neuron data */
void createSynapsePlot(vector<Synapse> synapse_output, ofstream &gpl, double t_max, string stimulus, double h_0, double theta_tag,
                       double theta_pro, double theta_p, double theta_d, int plh_cols)
{
	string x_max; // maximum x-value as string
	double x_scale = 100./96.; // value to enlarge x-axis
	string x_data; // string containing the column specifier for x-values
	string x_label = "Time "; // string containing the x-label

	if (t_max > 3600)
	{
		x_max = dtos(x_scale*t_max/60.);
		x_data = "($1/60)";
		x_label += "(min)";
	}
	else
	{
		x_max = dtos(x_scale*t_max, 1);
		x_data = "1";
		x_label += "(s)";
	}

	for (int i=0; i<synapse_output.size(); i++)
	{
		string spaces = "          ";
		string filename = string("synapse_plot_") + to_string(synapse_output[i].presyn_neuron) + string("_") + to_string(synapse_output[i].postsyn_neuron) + string(".gpl");
		ofstream f(filename);
		if (!f.is_open())
			return;

		f << "set term pdf enhanced font \"Sans, 18\"" << endl
		  << "set output \"" << dateStr("_synapse_") << synapse_output[i].presyn_neuron << "_" << synapse_output[i].postsyn_neuron << ".pdf\"" << endl
		  << endl
		  << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl
		  << endl
		  << "unset key" << endl // for better view
		  << "set label 1 at graph 0.04, graph 0.9 \"" << stimulus << "\" font \"Sans,22\" textcolor rgb rgb(255,123,0)" << endl
		  << "set ytics 1 nomirror" << endl
		  << endl
		  << "h0 = " << h_0 << endl
		  << "theta_pro = " << theta_pro << endl
		  << "theta_tag = " << theta_tag << endl
		  << "theta_p = " << theta_p << endl
		  << "theta_d = " << theta_d << endl
		  << endl
		  << "set xlabel \"" << x_label << "\"" << endl
		  << "set ylabel \"arb. units\"" << endl
		  << endl
		  << "set samples 1000" << endl
		  << "plot [0:" << x_max << "][0.0:3.2] \"" << dateStr("_data.txt") << "\" using " << x_data << ":(($" << plh_cols+i*3+2 <<
		     "+h0*$" << plh_cols+i*3+3 << ")/h0) title \"w(t)/h_0\" lw 2 lc rgb rgb(255,127,14) with lines, \\" << endl
		  << spaces << "\"" << dateStr("_data.txt") << "\" using " << x_data << ":" << plh_cols+i*3+4 << " title \"Ca\" lw 2 lc rgb rgb(200,200,150) with lines, \\" << endl
		  << spaces << "\"" << dateStr("_data.txt") << "\" using " << x_data << ":($" << plh_cols+i*3+2 << "/h0) title \"h(t)/h_0\" lw 2 lc rgb rgb(128,0,0) with lines, \\" << endl
		  << spaces << "\"" << dateStr("_data.txt") << "\" using " << x_data << ":($" << plh_cols+i*3+3 << "+1) title \"z(t)+1\" lw 2 lc rgb rgb(31,119,180) with lines, \\" << endl
		  << spaces << "1/0 title \" \" w dots lc rgb \"white\", \\" << endl
		  << spaces << "1/0 title \" \" w dots lc rgb \"white\", \\" << endl
		  << spaces << "1/0 title \" \" w dots lc rgb \"white\", \\" << endl
		  << spaces << "1/0 title \" \" w dots lc rgb \"white\", \\" << endl
		  << spaces << "1/0 title \" \" w dots lc rgb \"white\", \\" << endl
		  << spaces << "(h0+theta_pro)/h0 title \"{/Symbol q}_{pro}\" lc rgb rgb(0,128,0) lw 2 dt 2 with lines, \\" << endl
		  << spaces << "(h0+theta_tag)/h0 title \"{/Symbol q}_{tag}\" lc rgb rgb(255,0,0) lw 2 dt 2 with lines, \\" << endl
		  << spaces << "theta_p title \"{/Symbol q}_{p}\" lc rgb rgb(150,150,100) lw 2 dt 2 with lines, \\" << endl
		  << spaces << "theta_d title \"{/Symbol q}_{d}\" lc rgb rgb(150,150,150) lw 2 dt 2 with lines" << endl;

		// Call gnuplot script
		f.close();
		string gplc = concat("gnuplot ", filename);
		gpl << gplc << endl;
		system(gplc.c_str());
	}
}

/*** createExcNeuronPlot ***
 * Creates a plot that displays data specific to an excitatory neuron (voltage, protein amount) over time *
 * - neuron_output: a vector of the considered neurons *
 * - gpl: handle to a file containing the gnuplot calls *
 * - t_max: the duration of the time series in s */
void createExcNeuronPlot(vector<int> neuron_output, ofstream &gpl, double t_max)
{
	string x_max; // maximum x-value
	double x_scale = 100./96.; // value to enlarge x-axis
	string x_data; // string containing the column specifier for x-values
	string x_label = "Time "; // string containing the x-label

	if (t_max > 3600)
	{
		x_max = dtos(x_scale*t_max/60.);
		x_data = "($1/60)";
		x_label += "(min)";
	}
	else
	{
		x_max = dtos(x_scale*t_max, 1);
		x_data = "1";
		x_label += "(s)";
	}

	for (int i=0; i<neuron_output.size(); i++)
	{
		string spaces = "          ";
		string filename = string("exc_neuron_plot_") + to_string(neuron_output[i]) + string(".gpl");
		ofstream f(filename);
		if (!f.is_open())
			return;

		f << "set term pdf enhanced font \"Sans, 18\"" << endl
		  << "set output \"" << dateStr("_neuron_") << neuron_output[i] << ".pdf\"" << endl
		  << endl
		  << "unset key" << endl // for better view
		  << endl
		  << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl
		  << endl
		  << "set ytics nomirror" << endl
		  << "set y2tics tc rgb rgb(255,0,0)" << endl
		  << endl
		  << "set xlabel \"" << x_label << "\"" << endl
		  << "set ylabel \"p(t)\"" << endl
		  << "set y2label \"V(t) / mV\" tc rgb rgb(255,0,0)" << endl
		  << endl
		  << "set samples 1000" << endl
		  << "plot [0:" << x_max << "] \""
#if PROTEIN_POOLS == POOLS_C
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*3+4 << " axes x1y1 t \"Common protein\" lw 2 lc rgb rgb (0,128,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*3+2 << " axes x1y2 t \"Voltage\" lw 2 lc rgb rgb (255,0,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*3+3 << " axes x1y2 t \"Misc.\" lw 2 lc rgb rgb (255,128,0) with lines"
#elif PROTEIN_POOLS == POOLS_PD
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+4 << " axes x1y1 t \"LTP protein\" lw 2 lc rgb rgb (128,128,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+5 << " axes x1y1 t \"LTD protein\" lw 2 lc rgb rgb (0,128,128) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+2 << " axes x1y2 t \"Voltage\" lw 2 lc rgb rgb (255,0,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+3 << " axes x1y2 t \"Misc.\" lw 2 lc rgb rgb (255,128,0) with lines"

#elif PROTEIN_POOLS == POOLS_PCD
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+4 << " axes x1y1 t \"LTP protein\" lw 2 lc rgb rgb (128,128,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+5 << " axes x1y1 t \"Common protein\" lw 2 lc rgb rgb (0,128,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+6 << " axes x1y1 t \"LTD protein\" lw 2 lc rgb rgb (0,128,128) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+2 << " axes x1y2 t \"Voltage\" lw 2 lc rgb rgb (255,0,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << i*4+3 << " axes x1y2 t \"Misc.\" lw 2 lc rgb rgb (255,128,0) with lines"
#endif

//#if SPIKE_PLOTTING == NUMBER_AND_RASTER || SPIKE_PLOTTING == RASTER // plot spikes as dots - uses a lot of memory when showing the PDF
//		  << ", \\" << endl << spaces << "\""
//		  << dateStr("_spike_raster.txt") << "\" using 1:($2 == " << neuron_output[i] << " ? 0.5 : 1/0) axes x1y1 notitle with dots lc \"blue\""
//#endif
		  << endl;

		// Call gnuplot script
		f.close();
		string gplc = concat("gnuplot ", filename);
		gpl << gplc << endl;
		system(gplc.c_str());
	}
}

/*** createInhNeuronPlot ***
 * Creates a plot that displays the specific voltage of an inhibitory neuron over time *
 * - neuron_output: a vector of the considered neurons *
 * - gpl: handle to a file containing the gnuplot calls *
 * - t_max: the duration of the time series in s *
 * - plh_cols: number of columns to skip for they contain neuron data */
void createInhNeuronPlot(vector<int> neuron_output, ofstream &gpl, double t_max, int plh_cols)
{
	string x_max; // maximum x-value
	double x_scale = 100./96.; // value to enlarge x-axis
	string x_data; // string containing the column specifier for x-values
	string x_label = "Time "; // string containing the x-label

	if (t_max > 3600)
	{
		x_max = dtos(x_scale*t_max/60.);
		x_data = "($1/60)";
		x_label += "(min)";
	}
	else
	{
		x_max = dtos(x_scale*t_max, 1);
		x_data = "1";
		x_label += "(s)";
	}

	for (int i=0; i<neuron_output.size(); i++)
	{
		string spaces = "          ";
		string filename = string("inh_neuron_plot_") + to_string(neuron_output[i]) + string(".gpl");
		ofstream f(filename);
		if (!f.is_open())
			return;

		f << "set term pdf enhanced font \"Sans, 18\"" << endl
		  << "set output \"" << dateStr("_neuron_") << neuron_output[i] << ".pdf\"" << endl
		  << endl
		  << "unset key" << endl // for better view
		  << endl
		  << "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)" << endl
		  << endl
		  << "set xlabel \"" << x_label << "\"" << endl
		  << "set ylabel \"V(t) / mV\"" << endl
		  << endl
		  << "set samples 1000" << endl
		  << "plot [0:" << x_max << "] \"" << dateStr("_data.txt") << "\" using " << x_data << ":" << plh_cols+i*2+2 << " t \"Voltage\" lw 2 lc rgb rgb (255,0,0) with lines, \\" << endl << spaces << "\""
		  << dateStr("_data.txt") << "\" using " << x_data << ":" << plh_cols+i*2+3 
#if NEURON_MODEL == LIF
		  << " t \"Current\" "
#elif NEURON_MODEL == MAT2
		  << " t \"Threshold\" "
#endif
		  << "lw 2 lc rgb rgb (255,128,0) with lines"
#if SPIKE_PLOTTING == NUMBER_AND_RASTER || SPIKE_PLOTTING == RASTER
		  << ", \\" << endl << spaces << "\""
		  << dateStr("_spike_raster.txt") << "\" using 1:($2 == " << neuron_output[i] << " ? 0.5 : 1/0) axes x1y1 notitle with dots lc \"blue\""
#endif
		  << endl;

		// Call gnuplot script
		f.close();
		string gplc = concat("gnuplot ", filename);
		gpl << gplc << endl;
		system(gplc.c_str());
	}
}

/*** pyPlot ***
 * Runs a python function for plotting and adds the command to a script *
 * - cmd: the command to run the python function */
 
void pyPlot(string cmd)
{
	ofstream f ("pyplot", ofstream::out | ofstream::app);
	f << cmd << endl; // save command in script file
	f.close();
	
	system(cmd.c_str()); // execute command
}


/*** createNetworkPlotAveragedWeights ***
 * Creates a plot that displays firing rates and weights averaged over synapses; *
 * uses Python's Matplotlib and plotFunctions.py *
 * - t: the plot time in s *
 * - h_0: initial weight *
 * - Nl: number of neurons in one column/row *
 * - z_max: maximum late-phase coupling strength */
void createNetworkPlotAveragedWeights(double t, double h_0, int Nl, double z_max)
{
	string py_cmd = string("python3 -c \"import plotFunctions as pf; pf.plotAveragedWeights(") +
	                string("\'") + // filename
#if NETWORK_OUTPUT == COMPLETE
	                dateStr("_net_") + dtos(t,1) + string(".txt\', ") +
#elif NETWORK_OUTPUT == AVERAGED
	                dateStr("_net_av_") + dtos(t,1) + string(".txt\', ") +
#endif
	                dtos(h_0,NET_OUTPUT_PRECISION) + string(", ") + // h_0
	                string("-0.5, ") + // z_min
			dtos(z_max,1) + string(", ") + // z_max
	                dtos(Nl,0) + string(", ") + // Nl
#if NETWORK_OUTPUT == COMPLETE
	                string("False)\""); // already_averaged
#elif NETWORK_OUTPUT == AVERAGED
	                string("True)\""); // already_averaged
#endif

	pyPlot(py_cmd);
}

/*** createNetworkPlotWeights ***
 * Creates a plot that displays firing rates and weights; *
 * uses Python's Matplotlib and plotFunctions.py *
 * - t: the plot time in s *
 * - h_0: initial weight *
 * - Nl: number of neurons in one column/row *
 * - z_max: maximum late-phase coupling strength */
#if NETWORK_OUTPUT == COMPLETE
void createNetworkPlotWeights(double t, double h_0, int Nl, double z_max)
{
	string py_cmd = string("python3 -c \"import plotFunctions as pf; pf.plotWeights(") +
	                string("\'") + // filename
	                dateStr("_net_") + dtos(t,1) + string(".txt\', ") +
	                dtos(h_0,NET_OUTPUT_PRECISION) + string(", ") + // h_0
	                string("-0.5, ") + // z_min
			dtos(z_max,1) + string(", ") + // z_max
	                dtos(Nl,0) + string(")\""); // Nl
	pyPlot(py_cmd);
}
#endif

/*** plotMinSimResults ***
 * Creates a plot that displays crucial observables for one neuron and one synapse  *
 * - col_neur: the first column containing data of the targeted neuron (the membrane potential, the next two columns contain membrane current and, if applicable, protein amount) *
 * - col_syn: the first column containing data of the targeted synapse (the early-phase weight, the next two columns contain late-phase weight and calcium amount) *
 * - h_0: initial weight *
 * - theta_tag: tagging threshold *
 * - theta_pro: protein synthesis threshold *
 * - theta_p: potentiation threshold for Ca dynamics *
 * - theta_d: depression threshold for Ca dynamics *
 * - store_path: path to the resulting graphics file *
 */
void plotMinSimResults(int col_neur, int col_syn, double h_0, double theta_tag, double theta_pro, double theta_p, double theta_d, string store_path)
{
	string py_cmd = string("python3 -c \"import plotFunctions as pf; pf.plotMinOverview(\'")
	              + dateStr("_data.txt\'") + (",") 
	              + to_string(col_neur) + (",") + to_string(col_syn) + (",")
	              + dtos(h_0, OUTPUT_PRECISION) + (",") + dtos(theta_tag, OUTPUT_PRECISION) + (",")
	              + dtos(theta_pro, OUTPUT_PRECISION) + (",") 
	              + dtos(theta_p, OUTPUT_PRECISION) + (",") + dtos(theta_d, OUTPUT_PRECISION) + (",")
	              + ("'") + store_path + ("')\"");

	pyPlot(py_cmd);
}

