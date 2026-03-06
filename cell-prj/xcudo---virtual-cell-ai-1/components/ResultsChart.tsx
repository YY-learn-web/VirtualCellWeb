
import React, { useMemo } from 'react';
import { 
  BarChart, 
  Bar, 
  XAxis, 
  YAxis, 
  CartesianGrid, 
  Tooltip, 
  ResponsiveContainer,
  ReferenceLine,
  Cell,
  LineChart,
  Line,
  Legend
} from 'recharts';
import { GeneExpressionResult, SweepPoint } from '../types';

interface ResultsChartProps {
  data: GeneExpressionResult[];
  isDark: boolean;
}

export const ResultsChart: React.FC<ResultsChartProps> = ({ data, isDark }) => {
  
  // Sort data to show most upregulated and downregulated genes
  const sortedData = useMemo(() => {
    // Take top 20 upregulated and top 20 downregulated for clarity in demo
    const sorted = [...data].sort((a, b) => b.expressionLevel - a.expressionLevel);
    const top = sorted.slice(0, 20);
    const bottom = sorted.slice(sorted.length - 20, sorted.length);
    return [...top, ...bottom].sort((a, b) => b.expressionLevel - a.expressionLevel);
  }, [data]);

  // Colors aligned with Slate/Cyan theme
  const textColor = isDark ? '#94a3b8' : '#64748b'; // slate-400 : slate-500
  const gridColor = isDark ? '#1e293b' : '#e2e8f0'; // slate-800 : slate-200
  const tooltipBg = isDark ? '#0f172a' : '#ffffff'; // slate-900 : white
  const tooltipBorder = isDark ? '#334155' : '#e2e8f0'; // slate-700 : slate-200
  const tooltipText = isDark ? '#f8fafc' : '#0f172a'; // slate-50 : slate-900

  // Bar colors
  const posColor = isDark ? '#22d3ee' : '#0891b2'; // cyan-400 : cyan-600
  const negColor = isDark ? '#f472b6' : '#db2777'; // pink-400 : pink-600

  return (
    <div className="w-full h-[450px] bg-white dark:bg-slate-900 rounded-xl border border-slate-200 dark:border-slate-800 p-6 shadow-md shadow-slate-200/50 dark:shadow-none transition-all duration-300">
      <h3 className="text-sm font-mono text-slate-700 dark:text-cyan-400 mb-6 uppercase tracking-widest font-bold border-b border-slate-100 dark:border-slate-800 pb-2">
        Top Differential Genes (Prediction - Baseline)
      </h3>
      <ResponsiveContainer width="100%" height="85%">
        <BarChart
          data={sortedData}
          margin={{ top: 5, right: 30, left: 20, bottom: 5 }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke={gridColor} vertical={false} />
          <XAxis 
            dataKey="geneId" 
            tick={{ fill: textColor, fontSize: 10 }} 
            interval={0}
            angle={-45}
            textAnchor="end"
            height={60}
          />
          <YAxis 
            tick={{ fill: textColor, fontSize: 10 }}
            label={{ value: 'Expression Delta', angle: -90, position: 'insideLeft', fill: textColor, fontSize: 12 }} 
          />
          <Tooltip 
            contentStyle={{ backgroundColor: tooltipBg, borderColor: tooltipBorder, color: tooltipText, borderRadius: '8px', boxShadow: '0 10px 15px -3px rgb(0 0 0 / 0.1)' }}
            cursor={{ fill: isDark ? '#1e293b' : '#f1f5f9', opacity: 0.5 }}
          />
          <ReferenceLine y={0} stroke={isDark ? '#475569' : '#cbd5e1'} />
          <Bar dataKey="expressionLevel" radius={[4, 4, 0, 0]}>
            {sortedData.map((entry, index) => (
              <Cell 
                key={`cell-${index}`} 
                fill={entry.expressionLevel > 0 ? posColor : negColor} 
              />
            ))}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
};

interface SweepChartProps {
  data: SweepPoint[];
  xAxisLabel: string;
  isDark: boolean;
}

export const SweepLineChart: React.FC<SweepChartProps> = ({ data, xAxisLabel, isDark }) => {
  const textColor = isDark ? '#94a3b8' : '#64748b';
  const gridColor = isDark ? '#1e293b' : '#e2e8f0';
  const tooltipBg = isDark ? '#0f172a' : '#ffffff';
  const tooltipBorder = isDark ? '#334155' : '#e2e8f0';
  const tooltipText = isDark ? '#f8fafc' : '#0f172a';

  // Extract gene keys from the first data point
  const geneKeys = Object.keys(data[0]?.genes || {});
  
  // A palette of colors for the lines
  const colors = [
    '#22d3ee', // cyan
    '#a855f7', // purple
    '#f472b6', // pink
    '#22c55e', // green
    '#eab308'  // yellow
  ];

  // Flatten data for Recharts (move genes up to root level of object)
  const flatData = data.map(d => ({
    xValue: d.xValue,
    ...d.genes
  }));

  return (
     <div className="w-full h-[450px] bg-white dark:bg-slate-900 rounded-xl border border-slate-200 dark:border-slate-800 p-6 shadow-md shadow-slate-200/50 dark:shadow-none transition-all duration-300">
      <h3 className="text-sm font-mono text-slate-700 dark:text-cyan-400 mb-6 uppercase tracking-widest font-bold border-b border-slate-100 dark:border-slate-800 pb-2">
        Dose/Time Response Analysis (Top 5 Genes)
      </h3>
      <ResponsiveContainer width="100%" height="85%">
        <LineChart data={flatData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
          <CartesianGrid strokeDasharray="3 3" stroke={gridColor} vertical={false} />
          <XAxis 
            dataKey="xValue" 
            type="number"
            tick={{ fill: textColor, fontSize: 10 }}
            label={{ value: xAxisLabel, position: 'insideBottom', offset: -5, fill: textColor, fontSize: 12 }}
            domain={['dataMin', 'dataMax']}
          />
          <YAxis 
             tick={{ fill: textColor, fontSize: 10 }}
             label={{ value: 'Expression Delta', angle: -90, position: 'insideLeft', fill: textColor, fontSize: 12 }}
          />
          <Tooltip 
             contentStyle={{ backgroundColor: tooltipBg, borderColor: tooltipBorder, color: tooltipText, borderRadius: '8px' }}
          />
          <Legend wrapperStyle={{ paddingTop: '10px' }} />
          {geneKeys.map((gene, idx) => (
            <Line 
              key={gene}
              type="monotone" 
              dataKey={gene} 
              stroke={colors[idx % colors.length]} 
              strokeWidth={2}
              dot={{ r: 3, fill: colors[idx % colors.length] }}
            />
          ))}
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};
