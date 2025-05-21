
import React from 'react';

const Footer = () => {
  return (
    <footer className="bg-white border-t border-gray-100">
      <div className="container mx-auto py-12 px-4 sm:px-6 lg:px-8">
        <div className="grid grid-cols-1 md:grid-cols-4 gap-8">
          <div className="col-span-1 md:col-span-2">
            <div className="text-2xl font-bold text-qtl-blue mb-4">Attie Lab Resources</div>
            <p className="text-gray-600 mb-4 max-w-md">
              Advanced visualization and analysis tools for quantitative trait loci mapping in mice and other model organisms.
            </p>
            <div className="text-sm text-gray-500">
              © {new Date().getFullYear()} Attie Lab. All rights reserved.
            </div>
          </div>
          
          <div>
            <h3 className="text-gray-800 font-semibold mb-4">Resources</h3>
            <ul className="space-y-2 text-gray-600">
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Documentation</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Tutorials</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">API Reference</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Example Datasets</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Research Papers</a></li>
            </ul>
          </div>
          
          <div>
            <h3 className="text-gray-800 font-semibold mb-4">Lab</h3>
            <ul className="space-y-2 text-gray-600">
              <li><a href="#" className="hover:text-qtl-blue transition-colors">About Us</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Contact</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Privacy Policy</a></li>
              <li><a href="#" className="hover:text-qtl-blue transition-colors">Terms of Service</a></li>
            </ul>
          </div>
        </div>
      </div>
    </footer>
  );
};

export default Footer;
