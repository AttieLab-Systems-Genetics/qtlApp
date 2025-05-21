
import React from 'react';
import { Card, CardContent } from "@/components/ui/card";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";

const VisualizationDemo = () => {
  return (
    <section id="visualizations" className="py-20 bg-qtl-lightGray">
      <div className="container mx-auto px-4 sm:px-6 lg:px-8">
        <div className="text-center mb-16">
          <h2 className="text-3xl sm:text-4xl font-bold text-gray-900 mb-4">Interactive QTL Analysis</h2>
          <p className="text-xl text-gray-600 max-w-2xl mx-auto">
            Explore how our platform helps visualize and analyze QTL data with precision and clarity.
          </p>
        </div>

        <div className="max-w-5xl mx-auto">
          <Tabs defaultValue="lod-scores" className="w-full">
            <TabsList className="grid grid-cols-3 mb-8">
              <TabsTrigger value="lod-scores">LOD Score Plot</TabsTrigger>
              <TabsTrigger value="peaks-table">Peaks Table</TabsTrigger>
              <TabsTrigger value="cis-trans">Cis/Trans Analysis</TabsTrigger>
            </TabsList>
            
            <TabsContent value="lod-scores" className="mt-0">
              <Card>
                <CardContent className="p-6">
                  <div className="aspect-[16/9] bg-white rounded-lg overflow-hidden border border-gray-200">
                    <div className="h-full w-full p-4">
                      <div className="h-full flex flex-col">
                        <div className="flex justify-between items-center mb-2">
                          <div className="text-sm font-medium text-gray-500">Chromosome-wide LOD Score Plot</div>
                          <div className="flex gap-2">
                            <div className="px-2 py-1 text-xs bg-qtl-blue/10 text-qtl-blue rounded">Dataset: BXD</div>
                            <div className="px-2 py-1 text-xs bg-qtl-purple/10 text-qtl-purple rounded">Chr: 7</div>
                          </div>
                        </div>
                        <div className="flex-1 relative">
                          {/* Mock LOD score visualization inspired by screenshots */}
                          <div className="absolute inset-0">
                            <div className="h-full w-full relative">
                              {/* Grid background */}
                              <div className="absolute inset-0 grid grid-cols-12 grid-rows-6">
                                {Array.from({length: 72}).map((_, i) => (
                                  <div key={i} className="border-gray-100 border-[0.5px]"></div>
                                ))}
                              </div>
                              
                              {/* X-axis */}
                              <div className="absolute bottom-0 left-0 right-0 h-[1px] bg-gray-300"></div>
                              
                              {/* Y-axis */}
                              <div className="absolute top-0 bottom-0 left-0 w-[1px] bg-gray-300"></div>
                              
                              {/* Threshold line */}
                              <div className="absolute left-0 right-0 top-[30%] h-[1px] bg-red-400 border-t border-dashed border-red-500"></div>
                              <div className="absolute right-1 top-[30%] text-xs text-red-500 transform -translate-y-3">
                                Threshold (p=0.05)
                              </div>

                              {/* Y-axis labels */}
                              <div className="absolute left-1 top-0 text-xs text-gray-500">LOD</div>
                              <div className="absolute left-1 top-[75%] text-xs text-gray-500">2.0</div>
                              <div className="absolute left-1 top-[50%] text-xs text-gray-500">4.0</div>
                              <div className="absolute left-1 top-[25%] text-xs text-gray-500">6.0</div>

                              {/* X-axis labels */}
                              <div className="absolute bottom-1 left-[25%] text-xs text-gray-500">25 Mb</div>
                              <div className="absolute bottom-1 left-[50%] text-xs text-gray-500">50 Mb</div>
                              <div className="absolute bottom-1 left-[75%] text-xs text-gray-500">75 Mb</div>

                              {/* LOD curve with multiple peaks */}
                              <svg className="w-full h-full" preserveAspectRatio="none">
                                <path
                                  d="M0,85 C15,80 30,75 45,65 C60,55 75,45 90,35 C105,25 120,20 135,25 C150,35 165,55 180,75 C195,85 210,80 225,70 C240,60 255,50 270,30 C285,20 300,15 315,25 C330,40 345,60 360,70 C375,75 390,80 405,85 C420,88 435,85 450,80 C465,75 480,70 495,72 C510,75 525,80 540,82 C555,85 570,84 585,80"
                                  fill="none"
                                  stroke="#2563eb"
                                  strokeWidth="2.5"
                                  strokeLinecap="round"
                                  vectorEffect="non-scaling-stroke"
                                />
                              </svg>

                              {/* Peak markers */}
                              <div className="absolute left-[23%] bottom-0 top-[25%] w-[1px] bg-qtl-orange"></div>
                              <div className="absolute left-[23%] top-[25%] w-3 h-3 rounded-full bg-qtl-orange transform -translate-x-1.5 -translate-y-1.5"></div>
                              
                              <div className="absolute left-[45%] bottom-0 top-[15%] w-[1px] bg-qtl-orange"></div>
                              <div className="absolute left-[45%] top-[15%] w-4 h-4 rounded-full bg-qtl-orange transform -translate-x-2 -translate-y-2 ring-2 ring-qtl-orange/30"></div>
                              
                              <div className="absolute left-[68%] bottom-0 top-[20%] w-[1px] bg-qtl-orange"></div>
                              <div className="absolute left-[68%] top-[20%] w-3 h-3 rounded-full bg-qtl-orange transform -translate-x-1.5 -translate-y-1.5"></div>
                            </div>
                          </div>
                        </div>
                        <div className="mt-4 grid grid-cols-2 text-sm text-gray-600">
                          <div>
                            <span className="font-medium">Chromosome position</span>
                          </div>
                          <div className="text-right">
                            <span className="font-medium">Max LOD: 6.2</span>
                          </div>
                        </div>
                      </div>
                    </div>
                  </div>
                  <div className="mt-4 text-gray-600">
                    <p>Visualize LOD scores across chromosomes with our interactive plot. Easily identify significant peaks that exceed threshold values and zoom into regions of interest for detailed analysis.</p>
                  </div>
                </CardContent>
              </Card>
            </TabsContent>
            
            <TabsContent value="peaks-table" className="mt-0">
              <Card>
                <CardContent className="p-6">
                  <div className="aspect-[16/9] bg-white rounded-lg overflow-hidden border border-gray-200">
                    <div className="h-full w-full p-4">
                      <div className="h-full flex flex-col">
                        <div className="flex justify-between items-center mb-4">
                          <div className="text-sm font-medium text-gray-500">QTL Peaks Summary Table</div>
                          <div className="flex gap-2">
                            <div className="px-2 py-1 text-xs bg-qtl-blue/10 text-qtl-blue rounded">Dataset: BXD</div>
                            <div className="px-2 py-1 text-xs bg-qtl-purple/10 text-qtl-purple rounded">Threshold: 3.5</div>
                          </div>
                        </div>
                        <div className="flex-1 relative overflow-hidden">
                          {/* Mock peaks table */}
                          <div className="absolute inset-0 overflow-auto">
                            <table className="min-w-full divide-y divide-gray-200">
                              <thead className="bg-gray-50">
                                <tr>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Peak ID</th>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Chr</th>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Position (Mb)</th>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">LOD Score</th>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">p-value</th>
                                  <th scope="col" className="px-4 py-3 text-left text-xs font-medium text-gray-500 uppercase tracking-wider">Nearest Gene</th>
                                </tr>
                              </thead>
                              <tbody className="bg-white divide-y divide-gray-200">
                                <tr>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">QTL-001</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">7</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">23.4</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs font-medium text-qtl-blue">5.8</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">3.2e-06</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">Gata3</td>
                                </tr>
                                <tr>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">QTL-002</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">7</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">45.7</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs font-medium text-qtl-purple">6.2</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">8.5e-07</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">Sox9</td>
                                </tr>
                                <tr>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">QTL-003</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">7</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">68.2</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs font-medium text-qtl-blue">5.3</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">7.1e-06</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">Lef1</td>
                                </tr>
                                <tr>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">QTL-004</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">4</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">32.9</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs font-medium text-qtl-blue">4.9</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">1.3e-05</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">Foxp1</td>
                                </tr>
                                <tr>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">QTL-005</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">9</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">51.6</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs font-medium text-qtl-blue">4.2</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">4.8e-05</td>
                                  <td className="px-4 py-2 whitespace-nowrap text-xs">Runx2</td>
                                </tr>
                              </tbody>
                            </table>
                          </div>
                        </div>
                        <div className="mt-4 text-xs text-gray-500 text-right">
                          Showing 5 of 12 significant peaks
                        </div>
                      </div>
                    </div>
                  </div>
                  <div className="mt-4 text-gray-600">
                    <p>Our detailed peaks table summarizes all significant QTLs detected across chromosomes, providing positions, LOD scores, p-values, and nearest gene annotations for rapid biological interpretation.</p>
                  </div>
                </CardContent>
              </Card>
            </TabsContent>
            
            <TabsContent value="cis-trans" className="mt-0">
              <Card>
                <CardContent className="p-6">
                  <div className="aspect-[16/9] bg-white rounded-lg overflow-hidden border border-gray-200">
                    <div className="h-full w-full p-4">
                      <div className="h-full flex flex-col">
                        <div className="flex justify-between items-center mb-2">
                          <div className="text-sm font-medium text-gray-500">Cis/Trans eQTL Analysis</div>
                          <div className="flex gap-2">
                            <div className="px-2 py-1 text-xs bg-qtl-blue/10 text-qtl-blue rounded">Dataset: Liver eQTL</div>
                          </div>
                        </div>
                        <div className="flex-1 relative">
                          {/* Mock cis/trans plot visualization */}
                          <div className="absolute inset-0">
                            <div className="h-full w-full relative">
                              {/* Grid background */}
                              <div className="absolute inset-0 grid grid-cols-12 grid-rows-12">
                                {Array.from({length: 144}).map((_, i) => (
                                  <div key={i} className="border-gray-100 border-[0.5px]"></div>
                                ))}
                              </div>
                              
                              {/* Diagonal line for cis-eQTLs */}
                              <div className="absolute top-0 left-0 bottom-0 right-0">
                                <svg className="w-full h-full" preserveAspectRatio="none">
                                  <line 
                                    x1="0" 
                                    y1="100%" 
                                    x2="100%" 
                                    y2="0" 
                                    stroke="#94a3b8" 
                                    strokeWidth="1" 
                                    strokeDasharray="4" 
                                  />
                                </svg>
                              </div>

                              {/* Axis labels */}
                              <div className="absolute left-1 top-0 text-xs text-gray-500">Gene Position</div>
                              <div className="absolute left-0 bottom-1 text-xs text-gray-500">QTL Position</div>
                              
                              {/* Legend */}
                              <div className="absolute top-2 right-2 bg-white/80 p-1 rounded text-xs border border-gray-100">
                                <div className="flex items-center gap-1 mb-1">
                                  <div className="w-3 h-3 rounded-full bg-qtl-blue"></div>
                                  <span className="text-gray-600">Cis-eQTL</span>
                                </div>
                                <div className="flex items-center gap-1">
                                  <div className="w-3 h-3 rounded-full bg-qtl-orange"></div>
                                  <span className="text-gray-600">Trans-eQTL</span>
                                </div>
                              </div>

                              {/* Scatter plot dots - trans eQTLs */}
                              <div className="absolute top-[20%] left-[65%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[70%] left-[25%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[40%] left-[75%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[85%] left-[30%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[55%] left-[15%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[30%] left-[80%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              <div className="absolute top-[75%] left-[40%] w-2 h-2 bg-qtl-orange rounded-full"></div>
                              
                              {/* Scatter plot dots - cis eQTLs */}
                              <div className="absolute top-[80%] left-[15%] w-2.5 h-2.5 bg-qtl-blue rounded-full"></div>
                              <div className="absolute top-[60%] left-[35%] w-2.5 h-2.5 bg-qtl-blue rounded-full"></div>
                              <div className="absolute top-[45%] left-[50%] w-2.5 h-2.5 bg-qtl-blue rounded-full"></div>
                              <div className="absolute top-[30%] left-[65%] w-2.5 h-2.5 bg-qtl-blue rounded-full"></div>
                              <div className="absolute top-[15%] left-[80%] w-2.5 h-2.5 bg-qtl-blue rounded-full"></div>
                            </div>
                          </div>
                        </div>
                        <div className="mt-4 text-sm text-gray-600">
                          <div className="grid grid-cols-2 gap-4">
                            <div className="bg-qtl-blue/5 p-2 rounded">
                              <span className="font-medium text-qtl-blue">Cis-eQTLs:</span> 128
                            </div>
                            <div className="bg-qtl-orange/5 p-2 rounded">
                              <span className="font-medium text-qtl-orange">Trans-eQTLs:</span> 243
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>
                  </div>
                  <div className="mt-4 text-gray-600">
                    <p>Distinguish between cis and trans eQTL effects with our scatter plot visualization. Points near the diagonal represent cis-eQTLs where the genetic variant is close to the gene it affects, while distant points indicate trans-regulatory effects.</p>
                  </div>
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>
        </div>
      </div>
    </section>
  );
};

export default VisualizationDemo;
