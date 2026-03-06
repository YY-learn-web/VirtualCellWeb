import React from 'react';

interface InputProps extends React.InputHTMLAttributes<HTMLInputElement> {
  label: string;
  error?: string;
}

export const Input: React.FC<InputProps> = ({ label, error, className = '', ...props }) => {
  return (
    <div className="w-full mb-4">
      <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
        {label}
      </label>
      <input
        className={`w-full px-4 py-3 text-sm text-slate-900 dark:text-white bg-slate-50 dark:bg-slate-950/50 border rounded-lg focus:outline-none focus:ring-2 focus:ring-cyan-500/50 transition-all shadow-inner ${error ? 'border-red-500' : 'border-slate-200 dark:border-slate-700 focus:border-cyan-500'} ${className}`}
        {...props}
      />
      {error && <p className="mt-1 text-xs text-red-500 dark:text-red-400 font-mono">{error}</p>}
    </div>
  );
};

interface SliderProps extends React.InputHTMLAttributes<HTMLInputElement> {
  label: string;
  value: number;
  min: number;
  max: number;
  unit: string;
}

export const RangeSlider: React.FC<SliderProps> = ({ label, value, min, max, unit, ...props }) => {
  const percentage = ((value - min) / (max - min)) * 100;

  return (
    <div className="w-full mb-6">
      <div className="flex justify-between items-center mb-2">
        <label className="text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">{label}</label>
        <span className="text-sm font-bold text-slate-700 dark:text-white font-mono bg-slate-100 dark:bg-slate-800 px-2 py-1 rounded border border-slate-200 dark:border-slate-700">
          {value} <span className="text-slate-400 text-xs">{unit}</span>
        </span>
      </div>
      <div className="relative w-full h-2 bg-slate-200 dark:bg-slate-800 rounded-full">
        <div 
          className="absolute h-full bg-gradient-to-r from-cyan-600 to-cyan-400 rounded-full shadow-lg shadow-cyan-500/30" 
          style={{ width: `${percentage}%` }}
        />
        <input
          type="range"
          min={min}
          max={max}
          value={value}
          className="absolute w-full h-full opacity-0 cursor-pointer"
          {...props}
        />
      </div>
      <div className="flex justify-between mt-1 text-[10px] text-slate-400 dark:text-slate-500 font-mono">
        <span>{min}{unit}</span>
        <span>{max}{unit}</span>
      </div>
    </div>
  );
};