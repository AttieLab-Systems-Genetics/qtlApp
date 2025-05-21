
import React from 'react';
import { Dna, ArrowRight, ChartBar, Database, MessageSquare } from 'lucide-react';

const FeaturesSection = () => {
  const features = [
    {
      icon: <Dna className="h-10 w-10 text-qtl-blue" />,
      title: "Gene & Isoform QTL Analysis",
      description: "Explore gene expression QTLs and splicing variants with interactive visualizations and detailed LOD score plots."
    },
    {
      icon: <ChartBar className="h-10 w-10 text-qtl-blue" />,
      title: "Interactive LOD Score Plots",
      description: "Visualize QTL significance with interactive threshold adjustments and real-time identification of genetic loci."
    },
    {
      icon: <ArrowRight className="h-10 w-10 text-qtl-purple" />,
      title: "Liver & Clinical Trait Analysis",
      description: "Analyze liver lipids and clinical trait associations with specialized QTL mapping tools for metabolic phenotypes."
    },
    {
      icon: <ChartBar className="h-10 w-10 text-qtl-purple" />,
      title: "Phenotype Correlation Tools",
      description: "Discover relationships between clinical traits and genetic variants through correlation analysis and visualization."
    },
    {
      icon: <Database className="h-10 w-10 text-qtl-orange" />,
      title: "Founder Variant Database",
      description: "Access comprehensive founder strain genetic variant information through our searchable database portal."
    },
    {
      icon: <MessageSquare className="h-10 w-10 text-green-600" />,
      title: "Research Assistant Chatbot",
      description: "Get instant answers about Attie Lab research findings, experimental methods, and dataset information."
    }
  ];

  return (
    <section id="features" className="py-20 bg-white">
      <div className="container mx-auto px-4 sm:px-6 lg:px-8">
        <div className="text-center mb-16">
          <h2 className="text-3xl sm:text-4xl font-bold text-gray-900 mb-4">Attie Lab Research Tools</h2>
          <p className="text-xl text-gray-600 max-w-2xl mx-auto">
            Our applications provide researchers with comprehensive tools for investigating genetic factors in metabolic diseases.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
          {features.map((feature, index) => (
            <div 
              key={index} 
              className="bg-white p-6 rounded-xl border border-gray-100 shadow-sm hover:shadow-md transition-shadow duration-300"
            >
              <div className="bg-gray-50 p-3 rounded-lg inline-block mb-4">
                {feature.icon}
              </div>
              <h3 className="text-xl font-semibold text-gray-900 mb-3">{feature.title}</h3>
              <p className="text-gray-600">{feature.description}</p>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
};

export default FeaturesSection;
