
import React from 'react';
import { Button } from "@/components/ui/button";

const Header = () => {
  return (
    <header className="fixed w-full bg-white/90 z-50 border-b border-gray-100 shadow-sm">
      <div className="container mx-auto py-4 px-4 sm:px-6 lg:px-8">
        <div className="flex items-center justify-between">
          <div className="flex items-center">
            <span className="text-2xl font-bold text-qtl-blue">Attie Lab</span>
          </div>
          
          <nav className="hidden md:flex items-center space-x-8">
            <a href="#genes" className="text-gray-700 hover:text-qtl-blue transition-colors font-medium">Genes/Isoforms</a>
            <a href="#liver" className="text-gray-700 hover:text-qtl-purple transition-colors font-medium">Liver Lipids</a>
            <a href="#variants" className="text-gray-700 hover:text-qtl-orange transition-colors font-medium">Variant Portal</a>
            <a href="#chatbot" className="text-gray-700 hover:text-green-600 transition-colors font-medium">Chatbot</a>
          </nav>
          
          <div>
            <Button className="bg-qtl-blue hover:bg-blue-700 text-white">Get Started</Button>
          </div>
        </div>
      </div>
    </header>
  );
};

export default Header;
