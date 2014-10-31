package org.mutoss.gui;

import java.awt.Component;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.CancellationException;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.af.commons.widgets.ComponentTitledBorder;
import org.jdesktop.swingworker.SwingWorker;
import org.mutoss.config.ClassConfig;
import org.mutoss.config.Configuration;
import org.mutoss.gui.dialogs.TextFileViewer;

import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

public class AlgorithmPanel extends JPanel implements ActionListener, ChangeListener {
	
	JTextField jtTitle = new JTextField();
	JTextField jtReference = new JTextField();
	JTextField jtN1 = new JTextField("5000", 6);
	JTextField jtN2 = new JTextField("25", 6);
	JTextField jtRatio = new JTextField("1", 6);

	JButton ok = new JButton("Ready");
	JButton jbCompute = new JButton("Compute Design");
	public HTMLOutputPane jta;
	JLabel label = new JLabel();
	JPanel ntPanel = null;
	JPanel effPanel = null;
	List<JTextField> effV = new Vector<JTextField>();
	List<JTextField> nV = new Vector<JTextField>();
	CellConstraints cc = new CellConstraints();
	JRadioButton jbBalanceNothing = new JRadioButton("No balancing restrictions");	
	JRadioButton jbBalanceSequences = new JRadioButton("Balance treatments in regard to sequences");	 
	JRadioButton jbBalancePeriods = new JRadioButton("Balance treatments in regard to periods");
	CrossoverGUI gui;	
	JButton exportR = new JButton("Export to R");
	JButton showAlgoPerformance = new JButton("Search algorithm plot");
	JCheckBox useCatalogueDesigns = new JCheckBox("Use designs from catalogue as starting point");
	JComboBox jcbContrasts = new JComboBox(new String[] {"All pair comparisons (Tukey)", "Comparing treatment 1 to each of the others (Dunnett)"}); //, "User defined"} );
	String[] contrasts = new String[] {"Tukey", "Dunnett", "User defined"};
	//JComboBox jCBMixed = new JComboBox(new String[] {"Fixed subject effects model", "Random subject effects model"});
	JComboBox jcbCorrelation = new JComboBox(new String[] {"Independence", "Autoregressive Error", "Equicorrelated Error"}); //, "User defined"});
	String[] correlations = new String[] {"NULL", "autoregressive", "equicorrelated", "user defined"};
	JCheckBox fixedNumber = new JCheckBox("Specify exact number of treatment assignments:");
	JCheckBox fixedSubjectEffects = new JCheckBox("Include fixed subject effects in design matrix.");
	JLabel jlMixed;
	JLabel jlCor;
	JLabel jlVar;
	JTextField jtWithinSubjectRho;
	//JTabbedPane jTabAlgo = new jTabAlgo;
	//String udcm;
	JSplitPane pane;
	
	ClassConfig ac = new ClassConfig(Configuration.getInstance(), AlgorithmPanel.class);
	
	public AlgorithmPanel(CrossoverGUI gui) {
		this.gui = gui;
    	
    	jta = new HTMLOutputPane(gui);
		jta.setFont(new Font("Monospaced", Font.PLAIN, 12));
		//jta.setLineWrap(false);		
		//jta.setMargin(new Insets(4,4,4,4));
		jta.setEditable(false);		
		
		pane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(getLeftSidePanel()), getRightSidePanel());
		
		setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.BOTH;
		c.weightx=1; c.weighty=1;
		add(pane, c);
		
		loadDefaults();
		
		gui.spinnerT.addChangeListener(this);
		gui.spinnerP.addChangeListener(this);
		spinnerS.addChangeListener(this);
    	jcbCorrelation.addActionListener(this);
	}
	
	JSpinner spinnerS;
	JPanel lsPanel;

	public JPanel getRightSidePanel() {
		JPanel panel = new JPanel();
		
		String cols = "5dlu, pref, 5dlu, fill:min:grow, 5dlu";
        String rows = "5dlu, pref, 5dlu, fill:min:grow, 5dlu, pref, 5dlu";
        
        panel.setLayout(new FormLayout(cols, rows));
        CellConstraints cc = new CellConstraints();
		
		int row = 2;
    	
		panel.add(new JLabel("Created Design"), cc.xy(2, row));
		
		row += 2;
		
		panel.add(new JScrollPane(jta), cc.xyw(2, row, 3));
		jta.setEditable(true);
		
		row += 2;
		
		panel.add(exportR, cc.xy(2, row));
		exportR.setEnabled(false);
		exportR.addActionListener(this);		
		panel.add(showAlgoPerformance, cc.xy(4, row));
		showAlgoPerformance.setEnabled(false);
		showAlgoPerformance.addActionListener(this);
		
		return panel;
	}
	
	public void loadDefaults() {
		//JOptionPane.showMessageDialog(this, "Loading Settings");
		jcbCorrelation.setSelectedIndex(ac.getIntProperty("CVPattern", 0));
		jtWithinSubjectRho.setText(ac.getProperty("cpc", "0.5"));
		boolean mixed = jcbCorrelation.getSelectedIndex()==1 || jcbCorrelation.getSelectedIndex()==2;
		jtWithinSubjectRho.setEnabled(mixed);
		jlMixed.setEnabled(mixed);
		spinnerS.getModel().setValue(ac.getIntProperty("s", 4));		
		fixedSubjectEffects.setSelected(ac.getBoolProperty("fixedSubjectEffects", true));
		// effPanel is created more than one time. settings have to be loaded there.
		// TODO Could this result in inconsistencies?
		for (int i=0; i<nV.size(); i++) {
			nV.get(i).setText(ac.getProperty("nV"+i, ""));
		}
		jbBalanceNothing.setSelected(ac.getBoolProperty("jbBalanceNothing", true));
		jbBalanceSequences.setSelected(ac.getBoolProperty("jbBalanceSequences", false));
		jbBalancePeriods.setSelected(ac.getBoolProperty("jbBalancePeriods", false));
		jcbContrasts.setSelectedIndex(ac.getIntProperty("contrast", 0));
			
		useCatalogueDesigns.setSelected(ac.getBoolProperty("catalogue", false));
		jtN2.setText(ac.getProperty("nRuns", "20"));
		jtN1.setText(ac.getProperty("nSteps", "5000"));		
	}
    
	public void saveDefaults() {
		//JOptionPane.showMessageDialog(this, "Saving Settings");
		ac.setIntProperty("CVPattern", jcbCorrelation.getSelectedIndex());
		ac.setProperty("cpc", jtWithinSubjectRho.getText());
		ac.setIntProperty("s", Integer.parseInt(spinnerS.getModel().getValue().toString()));
		ac.setBoolProperty("fixedNumber", fixedNumber.isSelected());
		ac.setBoolProperty("fixedSubjectEffects", fixedSubjectEffects.isSelected());		
		for (int i=0; i<nV.size(); i++) {
			ac.setProperty("nV"+i, nV.get(i).getText());			
		}
		ac.setBoolProperty("jbBalanceNothing", jbBalanceNothing.isSelected());
		ac.setBoolProperty("jbBalanceSequences", jbBalanceSequences.isSelected());
		ac.setBoolProperty("jbBalancePeriods", jbBalancePeriods.isSelected());
		ac.setIntProperty("contrast", jcbContrasts.getSelectedIndex());
		for (int i=0; i<effV.size(); i++) {
			ac.setProperty("effV"+i, effV.get(i).getText());			
		}
		ac.setBoolProperty("catalogue", useCatalogueDesigns.isSelected());		
		ac.setProperty("nRuns", jtN2.getText());
		ac.setProperty("nSteps", jtN1.getText());
	}
	
	public JPanel getLeftSidePanel() {
		lsPanel = new JPanel();
		String cols = "5dlu, pref, 5dlu, fill:min:grow, 5dlu";
        String rows = "5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu";
        
        lsPanel.setLayout(new FormLayout(cols, rows));
		
		int row = 2;
    	
		/*

        lsPanel.add(jCBMixed, cc.xy(4, row));
        jCBMixed.setSelectedIndex(0);
        jCBMixed.addActionListener(this);
        //jCBMixed.setEnabled(false);
        
        row+=2;  
        
        jlVar = new JLabel("Random subject Var divided by Var of ε"); 
        jlVar.setEnabled(false);
        jtRatio.setEnabled(false);
        lsPanel.add(jlVar, cc.xy(2, row));
		lsPanel.add(jtRatio, cc.xy(4, row));
        
		row+=2;  
		
		*/
		
        jlCor = new JLabel("Covariance pattern");
        //jlCor.setEnabled(false);
        //jcbCorrelation.setEnabled(false);
		lsPanel.add(jlCor, cc.xy(2, row));
		lsPanel.add(jcbCorrelation, cc.xy(4, row));
		

        row+=2;
        
        jlMixed = new JLabel("Covariance pattern coefficient");
        jlMixed.setEnabled(false);
        jtWithinSubjectRho = new JTextField("0.5");
        jtWithinSubjectRho.setEnabled(false);
		lsPanel.add(jlMixed, cc.xy(2, row));
		lsPanel.add(jtWithinSubjectRho, cc.xy(4, row));
		
		row+=2;
		
		lsPanel.add(fixedSubjectEffects, cc.xyw(2, row, 3));		
		
        row+=2;  
        
    	spinnerS = new JSpinner(new SpinnerNumberModel(5, 1, 100, 1)); 
    	        
        lsPanel.add(new JLabel("Number of sequences:"), cc.xy(2, row));
        lsPanel.add(spinnerS, cc.xy(4, row));
        
        row+=2;
        
        rowN = row; 
        createTreatmentNumberPanel(true);
        
        row+=2;  
        
        ButtonGroup group = new ButtonGroup();
        group.add(jbBalanceNothing);
        group.add(jbBalanceSequences);
        jbBalanceSequences.setToolTipText("May decrease efficiency");
        group.add(jbBalancePeriods);
        jbBalancePeriods.setToolTipText("May decrease efficiency");
        jbBalanceNothing.setSelected(true);
              
        lsPanel.add(jbBalanceNothing, cc.xyw(2, row, 3));
        row+=2;
        lsPanel.add(jbBalanceSequences, cc.xyw(2, row, 3));
        row+=2;
        lsPanel.add(jbBalancePeriods, cc.xyw(2, row, 3));        
        row+=2;
        
        lsPanel.add(new JLabel("Contrasts:"), cc.xy(2, row));
        
        row+=2;
        
        lsPanel.add(jcbContrasts, cc.xyw(2, row, 3));

        row+=2;    
        
        rowEff = row;        
        createEffPanel();
        
        row+=2;        
        
        useCatalogueDesigns.setSelected(true);
        lsPanel.add(useCatalogueDesigns, cc.xyw(2, row, 3));
        
        row+=2;     
        
        lsPanel.add(new JLabel("Number of search runs:"), cc.xy(2, row));
        lsPanel.add(jtN2, cc.xy(4, row));
        
        row+=2;   
        
        lsPanel.add(new JLabel("Number of steps per run:"), cc.xy(2, row));
        lsPanel.add(jtN1, cc.xy(4, row));
        
        row+=2;        
        
        lsPanel.add(jbCompute, cc.xyw(2, row, 3));
        jbCompute.addActionListener(this);
        
        return lsPanel;
	}
	
	private int rowN, rowEff;
	
	public void createEffPanel() {		
		if (effPanel!=null) {
			lsPanel.remove(effPanel);
		}

		effPanel = new JPanel();
		GridBagConstraints effWeights = new GridBagConstraints();
		effWeights.fill = GridBagConstraints.BOTH;	
		effWeights.gridx=0; effWeights.gridy=0;
		effWeights.gridwidth = 1; effWeights.gridheight = 1;
		effWeights.ipadx=5; effWeights.ipady=5;
		effWeights.weightx=1; effWeights.weighty=1;
		effPanel.setBorder(BorderFactory.createTitledBorder("Weights:"));
		
		//effPanel.setEnabled(false);
		effPanel.setLayout(new GridBagLayout());
		
		List<String> labels = new Vector<String>();
        
		int s = gui.jCBmodel.getSelectedIndex();
		for (String p : CrossoverGUI.parameters[s]) {        	
        	labels.add(p);        	
        }  

		effV.clear();
		for (int i=0;i<labels.size();i++) {        		
			effV.add(new JTextField("1.0", 6));
			effPanel.add(new JLabel(labels.get(i)), effWeights);
			effWeights.gridx++;
			effPanel.add(effV.get(i), effWeights);	
			if ((i+1)%3!=0) {
				effWeights.gridx++;
			} else {
				effWeights.gridx=0;effWeights.gridy++;
			}
		}	
		if (labels.size()==0) {
			effPanel.add(new JLabel("Not applicable."), effWeights);
		}
		for (int i=0; i<effV.size(); i++) {
			effV.get(i).setText(ac.getProperty("effV"+i, i==0?"1":"0"));
		}	
		//for (Component c : effPanel.getComponents()) c.setEnabled(false);
		lsPanel.add(effPanel, cc.xyw(2, rowEff, 3));
		lsPanel.revalidate();
		lsPanel.repaint();
	}
	
	private String getEffFactors() {
		if (effV.size()==0) return "1";
		String s="c(";
		for (JTextField jt : effV) {
			s = s + jt.getText()+",";
		}
		return s.substring(0, s.length()-1)+")";
	}
	
	private void createTreatmentNumberPanel(boolean firstCall) {
		if (ntPanel!=null) {
			lsPanel.remove(ntPanel);
		}

		ntPanel = new JPanel();
		GridBagConstraints cbcWeights = new GridBagConstraints();
		cbcWeights.fill = GridBagConstraints.BOTH;	
		cbcWeights.gridx=0; cbcWeights.gridy=0;
		cbcWeights.gridwidth = 1; cbcWeights.gridheight = 1;
		//cbcWeights.ipadx=5; cbcWeights.ipady=5;
		cbcWeights.weightx=1; cbcWeights.weighty=1;
		cbcWeights.insets = new Insets(2, 2, 2, 2);
		if (firstCall) fixedNumber.addActionListener(this);
		ntPanel.setBorder(new ComponentTitledBorder(fixedNumber, ntPanel, BorderFactory.createTitledBorder("Weights:")));
		//ntPanel.setBorder(BorderFactory.createTitledBorder("Number of treatment assignments"));

		ntPanel.setLayout(new GridBagLayout());
		
		List<String> labels = new Vector<String>();
        
		int v = Integer.parseInt(gui.spinnerT.getModel().getValue().toString());
		int s = Integer.parseInt(spinnerS.getModel().getValue().toString());
		int p = Integer.parseInt(gui.spinnerP.getModel().getValue().toString());
		
		for (int i=1; i<=v; i++) {        	
        	labels.add(" Treatment "+i+":");        	        	
        }      
        
		nV.clear();
		for (int i=0;i<labels.size();i++) {        		
			nV.add(new JTextField(""+((s*p)/v+((i<(s*p)%v)?1:0)), 3));
			ntPanel.add(new JLabel(labels.get(i)), cbcWeights);
			cbcWeights.gridx++;
			ntPanel.add(nV.get(i), cbcWeights);	
			if ((i+1)%3!=0) {
				cbcWeights.gridx++;
			} else {
				cbcWeights.gridx=0;cbcWeights.gridy++;
			}
		}		
		boolean fixed = ac.getBoolProperty("fixedNumber", false);
		fixedNumber.setSelected(fixed);
		ntPanel.setEnabled(fixed);
		for (Component c : ntPanel.getComponents()) c.setEnabled(fixed);
		lsPanel.add(ntPanel, cc.xyw(2, rowN, 3));
		lsPanel.revalidate();
		lsPanel.repaint();
	}
	
	private String getVRep() {
		if (!fixedNumber.isSelected()) return "";
		String vrep = ", v.rep=c(";
		for (int i=0; i<nV.size(); i++) {   
			vrep += nV.get(i).getText();
			if (i!=nV.size()-1) vrep += ",";
		}
		return vrep+")";
	}
	
	private String getContrast() {
		return "";
	}
	
	public void addActionListener(ActionListener al) {
		ok.addActionListener(al);
	}
	
	public Design getDesign() {
		return null;
	}
	
	String command = "";
	String models = "";	

	public void actionPerformed(ActionEvent e) {
		saveDefaults();
		//System.out.println(""+e.getActionCommand()+", "+e.toString());
		if (e.getSource() == jbCompute) {
			gui.glassPane.start();
			
			command = "searchCrossOverDesign(s="+spinnerS.getModel().getValue().toString()
					+", "+gui.getParameters()
					+", model=\""+gui.jCBmodel.getSelectedItem()+"\""
					+", eff.factor="+getEffFactors()
					+", contrast=\""+contrasts[jcbContrasts.getSelectedIndex()]+"\""
					+getVRep()
					+", balance.s="+(jbBalanceSequences.isSelected()?"TRUE":"FALSE")
					+", balance.p="+(jbBalancePeriods.isSelected()?"TRUE":"FALSE")
					+getCorrelation()
					+(gui.jCBmodel.getSelectedIndex()==gui.PLACEBOMODEL?", model.param=list(placebos="+gui.jtfParam.getText()+")":"")
					+(gui.jCBmodel.getSelectedIndex()==gui.PROPORTIONALMODEL?", model.param=list(ppp="+gui.jtfParam.getText()+")":"")
					+(useCatalogueDesigns.isSelected()?", start.designs=\"catalog\"":"")	
					+", random.subject="+(fixedSubjectEffects.isSelected()?"FALSE":"TRUE")
					+", n=c("+jtN1.getText()+","+jtN2.getText()+")"
					+", verbose=FALSE)";
			
			//startTesting();		
			SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {
				//RList result;
				//String table;

				@Override
				protected Void doInBackground() throws Exception {					
					RControl.getR().eval(".COresult <- "+command);
					//RControl.getR().eval(".COresult <- getDesign(.COresult)");
					//table = RControl.getR().eval("Crossover:::getTable(getDesign(.COresult))").asRChar().getData()[0];
					return null;
				}

				protected final void done() {
					try {
						Design design = new Design("Search Result n=("+jtN1.getText()+","+jtN2.getText()+"), model="+gui.jCBmodel.getSelectedIndex()+", "+(new SimpleDateFormat("yyyy-MM-dd kk:mm").format(new Date())), ".COresult", command);
						gui.dac.addSearchResult(design);
						get();
						exportR.setEnabled(true);						
						showAlgoPerformance.setEnabled(true);
						jta.clear();
						jta.showDesign(design);
						jta.appendParagraph("Random seed: TODO");
						jta.appendParagraph("R Code: <pre>"+command+"</pre>");
						jta.textArea.setCaretPosition(0);
					} catch (CancellationException e) {
						// Will be perhaps used in the future.
					} catch (Throwable e) {
						String message = e.getMessage();
						//System.out.println("\""+message+"\"");
						if (message.equals("Error: \n")) message = "Empty message (most likely an error in the C++ code - please look at the R console for further output)\n\n";
						JOptionPane.showMessageDialog(gui, "R call produced an error:\n\n"+message+"\nWe will open a window with R code to reproduce this error for investigation.", "Error in R Call", JOptionPane.ERROR_MESSAGE);
						String traceback = RControl.getR().eval("paste(unlist(traceback()),collapse=\"\\n\")").asRChar().getData()[0];
						JDialog d = new JDialog(gui, "R Error", true);
						d.add( new TextFileViewer(gui, "R Objects", "The following R code produced the following error:\n\n" +message+"\n\nTraceback:\n\n"+traceback+"\n\n"+
										command, true) );
						d.pack();
						d.setSize(800, 600);
						d.setVisible(true);
						e.printStackTrace();
					} finally {
						gui.glassPane.stop();
					}
				}				
			};
			worker.execute();
			
			/*
			command = "s="+spinnerS.getModel().getValue().toString()
					+", "+gui.getParameters()
					+", model=\""+gui.jCBmodel.getSelectedItem()+"\""
					//+", eff.factor="+1					
					+(fixedNumber.isSelected()?", v.rep="+getVRep():"")
					+", balance.s="+(jbBalanceSequences.isSelected()?"TRUE":"FALSE")
					+", balance.p="+(jbBalancePeriods.isSelected()?"TRUE":"FALSE")
					+(gui.jCBmodel.getSelectedIndex()==4?", model.param=list(placebos="+gui.jtfParam.getText()+")":"")
					+(gui.jCBmodel.getSelectedIndex()==7?", model.param=list(ppp="+gui.jtfParam.getText()+")":"")
					+", verbose=FALSE"
					+", n=c("+jtN1.getText()+", 5)"
					;
			models = (useCatalogueDesigns.isSelected()?", start.designs=\"catalog\"":"");
			
			InfiniteRunningDialog ird = new InfiniteRunningDialog(gui, command, models);
			*/
		} else if (e.getSource()==showAlgoPerformance) {			
			RControl.getR().eval("JavaGD(\"Search algorithm performance\")");
			//RControl.getR().eval("png(filename=\""+getTmpFile()+"\")");
			RControl.getR().eval("print(plot(.COresult))");
			//RControl.getR().eval("dev.off()");
		} else if (e.getSource()==exportR) {
			VariableNameDialog vd = new VariableNameDialog(gui, "design");			
			RControl.getR().eval(vd.getName()+" <- getDesign(.COresult)");
			/*} else if (e.getSource()==jCBMixed) {
			boolean mixed = jCBMixed.getSelectedIndex()==1;
	        jlVar.setEnabled(mixed);
	        jtRatio.setEnabled(mixed);*/
		} else if (e.getSource()==jcbCorrelation) {
			boolean mixed = jcbCorrelation.getSelectedIndex()==1 || jcbCorrelation.getSelectedIndex()==2;
			jtWithinSubjectRho.setEnabled(mixed);
			jlMixed.setEnabled(mixed);
		} else if (e.getSource()==fixedNumber) {
			boolean fixed = fixedNumber.isSelected();
			ntPanel.setEnabled(fixed);
			for (Component c : ntPanel.getComponents()) c.setEnabled(fixed);
		}
	}

	private String getCorrelation() {
		if (jcbCorrelation.getSelectedIndex()==0) return "";
		//if (jcbCorrelation.getSelectedIndex()==3) return ", correlation="+udcm;		
		return ", correlation=\""+correlations[jcbCorrelation.getSelectedIndex()]+"\", rho="+jtWithinSubjectRho.getText();
	}


	private String getTmpFile() {
		String file = System.getProperty("java.io.tmpdir")+"/searchplot"+(new Date()).getTime()+".png";
		System.out.println(file);
		return file;
	}


	public void stateChanged(ChangeEvent e) {
		createTreatmentNumberPanel(false);
		int s1 = Integer.parseInt(gui.spinnerS1.getModel().getValue().toString());
		int s2 = Integer.parseInt(gui.spinnerS2.getModel().getValue().toString());
		int s = Integer.parseInt(spinnerS.getModel().getValue().toString());
		if (s<s1) {
			gui.spinnerS1.getModel().setValue(s);
		} else if (s>s2) {
			gui.spinnerS2.getModel().setValue(s);
		}
	}


	public void searchResultReady(String table) {
		exportR.setEnabled(true);						
		showAlgoPerformance.setEnabled(true);
		jta.clear();
		//jta.appendHTML(table);
		String command2 = "paste(capture.output(general.carryover(.COresult)),collapse=\"\\n\")";
		jta.appendParagraph("<pre>"+RControl.getR().eval(command2).asRChar().getData()[0]+"</pre>");		
		jta.appendParagraph("Random seed: TODO");
		jta.appendParagraph("R Code: <pre>"+command+"</pre>");
		//jta.setCaretPosition(0);		
	}

}
