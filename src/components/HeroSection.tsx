
import React from 'react';
import { Button } from "@/components/ui/button";
import { Dna, ArrowRight, Database, MessageSquare } from 'lucide-react';
import { Card, CardContent } from "@/components/ui/card";

const HeroSection = () => {
  return (
    <section className="pt-24 pb-20 bg-gradient-to-br from-white via-slate-50 to-gray-100 relative overflow-hidden">
      {/* Simplified background */}
      <div className="absolute inset-0 z-0 opacity-10">
        <div className="absolute top-20 left-10 w-72 h-72 rounded-full bg-qtl-blue/20 blur-3xl"></div>
        <div className="absolute bottom-20 right-10 w-80 h-80 rounded-full bg-qtl-purple/20 blur-3xl"></div>
      </div>
      
      <div className="container mx-auto px-4 sm:px-6 lg:px-8 relative z-10">
        <div className="text-center mb-16">
          <h1 className="text-4xl sm:text-6xl font-bold leading-tight mb-8 text-gray-800">
            Welcome to <span className="bg-gradient-to-r from-qtl-blue to-qtl-purple bg-clip-text text-transparent">Attie Lab</span> Resources
          </h1>
          <p className="text-xl text-gray-700 max-w-3xl mx-auto">
            Advanced tools for genetics research, QTL mapping, and data analysis to accelerate discoveries in metabolic diseases.
          </p>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-8 max-w-5xl mx-auto">
          <Card className="border border-blue-100 shadow-md hover:shadow-lg transition-shadow duration-300">
            <CardContent className="p-8">
              <div className="bg-qtl-blue/10 p-4 rounded-full mb-6 inline-flex">
                <Dna size={40} className="text-qtl-blue" />
              </div>
              <h3 className="text-2xl font-semibold mb-4 text-gray-800">Liver Gene & Liver Isoforms QTL</h3>
              <p className="text-gray-600 mb-6 text-lg">Analyze liver gene expression and splicing variants with advanced QTL mapping tools.</p>
              <Button className="mt-auto bg-qtl-blue hover:bg-blue-700 text-white w-full">
                <Dna className="mr-2" /> Explore Liver Gene QTLs
              </Button>
            </CardContent>
          </Card>
          
          <Card className="border border-purple-100 shadow-md hover:shadow-lg transition-shadow duration-300">
            <CardContent className="p-8">
              <div className="bg-qtl-purple/10 p-4 rounded-full mb-6 inline-flex">
                <ArrowRight size={40} className="text-qtl-purple" />
              </div>
              <h3 className="text-2xl font-semibold mb-4 text-gray-800">Liver Lipids & Clinical Traits</h3>
              <p className="text-gray-600 mb-6 text-lg">Investigate liver lipids and clinical traits with specialized phenotype QTL analysis.</p>
              <Button className="mt-auto bg-qtl-purple hover:bg-purple-700 text-white w-full">
                <ArrowRight className="mr-2" /> Explore Liver QTLs
              </Button>
            </CardContent>
          </Card>
          
          <Card className="border border-orange-100 shadow-md hover:shadow-lg transition-shadow duration-300">
            <CardContent className="p-8">
              <div className="bg-qtl-orange/10 p-4 rounded-full mb-6 inline-flex">
                <Database size={40} className="text-qtl-orange" />
              </div>
              <h3 className="text-2xl font-semibold mb-4 text-gray-800">Founder Variant Portal</h3>
              <p className="text-gray-600 mb-6 text-lg">Access comprehensive founder strain variant data through our interactive portal.</p>
              <Button className="mt-auto bg-qtl-orange hover:bg-orange-700 text-white w-full">
                <Database className="mr-2" /> Browse Variants
              </Button>
            </CardContent>
          </Card>
          
          <Card className="border border-green-100 shadow-md hover:shadow-lg transition-shadow duration-300">
            <CardContent className="p-8">
              <div className="bg-green-100 p-4 rounded-full mb-6 inline-flex">
                <MessageSquare size={40} className="text-green-600" />
              </div>
              <h3 className="text-2xl font-semibold mb-4 text-gray-800">Attie Lab Chatbot</h3>
              <p className="text-gray-600 mb-6 text-lg">Get instant answers about our datasets, experimental methods, and research findings.</p>
              <Button className="mt-auto bg-green-600 hover:bg-green-700 text-white w-full">
                <MessageSquare className="mr-2" /> Chat with Attie
              </Button>
            </CardContent>
          </Card>
        </div>
      </div>
    </section>
  );
};

export default HeroSection;
