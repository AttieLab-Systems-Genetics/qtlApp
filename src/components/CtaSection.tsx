
import React from 'react';
import { Button } from "@/components/ui/button";

const CtaSection = () => {
  return (
    <section id="about" className="py-20 bg-gradient-to-br from-qtl-blue/10 to-qtl-purple/10">
      <div className="container mx-auto px-4 sm:px-6 lg:px-8">
        <div className="max-w-4xl mx-auto text-center">
          <h2 className="text-3xl sm:text-4xl font-bold text-gray-900 mb-6">Ready to Elevate Your Genetic Research?</h2>
          <p className="text-xl text-gray-600 mb-10 max-w-2xl mx-auto">
            Join researchers worldwide who are using our QTL App to accelerate discoveries in quantitative trait analysis and genetic mapping.
          </p>
          
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8 mb-12">
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-100">
              <div className="text-3xl font-bold text-qtl-blue mb-2">30+</div>
              <div className="text-gray-600">Research Institutions</div>
            </div>
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-100">
              <div className="text-3xl font-bold text-qtl-purple mb-2">250+</div>
              <div className="text-gray-600">Published Studies</div>
            </div>
            <div className="bg-white p-6 rounded-xl shadow-sm border border-gray-100">
              <div className="text-3xl font-bold text-qtl-orange mb-2">10k+</div>
              <div className="text-gray-600">Analyses Run</div>
            </div>
          </div>
          
          <div className="space-y-6">
            <div className="flex flex-col sm:flex-row justify-center gap-4">
              <Button className="bg-qtl-blue hover:bg-blue-700 text-white px-8 py-6 text-lg">
                Start Free Trial
              </Button>
              <Button variant="outline" className="border-qtl-purple text-qtl-purple hover:bg-qtl-purple/5 px-8 py-6 text-lg">
                Schedule a Demo
              </Button>
            </div>
            <p className="text-sm text-gray-500">
              No credit card required. Free 14-day trial with full access to all features.
            </p>
          </div>
        </div>
      </div>
    </section>
  );
};

export default CtaSection;
