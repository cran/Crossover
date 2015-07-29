package org.mutoss.gui;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.lang.reflect.Method;
import java.text.DecimalFormat;
import java.util.Locale;
import java.util.Properties;

import javax.swing.ImageIcon;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.af.commons.errorhandling.DefaultExceptionHandler;
import org.af.commons.errorhandling.ErrorHandler;
import org.af.commons.logging.ApplicationLog;
import org.af.commons.logging.LoggingSystem;
import org.af.commons.tools.OSTools;
import org.af.commons.widgets.InfiniteProgressPanel;
import org.af.commons.widgets.InfiniteProgressPanel.AbortListener;
import org.af.commons.widgets.WidgetFactory;
import org.apache.log4j.PropertyConfigurator;
import org.jdesktop.swingworker.SwingWorker;
import org.mutoss.config.ClassConfig;
import org.mutoss.config.Configuration;
import org.mutoss.gui.archive.DesignArchiveControl;
import org.mutoss.gui.dialogs.ConfirmArchivePath;
import org.mutoss.gui.dialogs.ErrorDialogChooseLevel;
import org.mutoss.gui.dialogs.ErrorDialogSGTK;
import org.mutoss.gui.dialogs.TellAboutOnlineUpate;

import com.jgoodies.forms.layout.CellConstraints;
import com.jgoodies.forms.layout.FormLayout;

public class CrossoverGUI extends JFrame implements WindowListener, ActionListener, ChangeListener, ListSelectionListener, AbortListener, DocumentListener {
	 static { // Static initializer block to set the Locale
		 Locale.setDefault(Locale.ENGLISH);
		 JComponent.setDefaultLocale(Locale.ENGLISH);			 
	 };
	 
	Configuration conf = Configuration.getInstance();
	public InfiniteProgressPanel glassPane;
	DesignArchiveControl dac;
	
	public static DecimalFormat df = new DecimalFormat("#.####");
	
	public static final String[] models = new String[] {
			"Standard additive model",
			"Self-adjacency model",
			"Proportionality model",
			"Placebo model",
			"No carry-over into self model",
			"Treatment decay model",
			"Full set of interactions",
			"Second-order carry-over effects",
			"No carry-over effects"
	};
	
	public static final String[][] parameters = new String[][] {
		{"Treatment", "Carry-over"},
		{"Treatment", "Carry-over", "Carry-over (into-self)"},
		{},
		{"Treatment", "Carry-over"},
		{"Treatment", "Carry-over"},
		{"Treatment", "Carry-over"},		
		{"Treatment", "Carry-over", "Interaction"},
		{"Treatment", "Carry-over (1st order)", "Carry-over (2nd order)"},
		{}
	};
	
	public static final int PLACEBOMODEL = 3;
	public static final int PROPORTIONALMODEL = 2;
	 
	/**
	 * Constructor - you will not need to change things here
	 */
	public CrossoverGUI() {		
		super("CROSS-OVER DESIGN SEARCH TOOL");
		//setIconImage((new ImageIcon(getClass().getResource("/org/mutoss/gui/graph/images/rjavaicon64.png"))).getImage());
		
		if (Configuration.getInstance().getGeneralConfig().failSafeMode() 
				&& JOptionPane.showConfirmDialog(this, "Reset settings?", "Reset settings?", JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION) {
			Configuration.getInstance().clearConfiguration();
		}
		
		if (!Configuration.getInstance().getGeneralConfig().failSafeMode() 
				|| JOptionPane.showConfirmDialog(this, "Start Logging and Error Handler?", "Logging and error handler", JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION) {

			Properties log4jProperties = new Properties();
			log4jProperties.setProperty("log4j.rootLogger", "DEBUG, READABLE_FILE");
			log4jProperties.setProperty("log4j.appender.READABLE_FILE", "org.apache.log4j.FileAppender");
			log4jProperties.setProperty("log4j.appender.READABLE_FILE.layout", "org.apache.log4j.PatternLayout");
			log4jProperties.setProperty("log4j.appender.READABLE_FILE.layout.conversionPattern", "%-5p - %m%n       [%t] (%c:%M at %F:%L)%n");
			log4jProperties.setProperty("log4j.appender.READABLE_FILE.File", "/tmp/mutoss.log");
			log4jProperties.setProperty("log4j.appender.READABLE_FILE.append", "false");
			PropertyConfigurator.configure(log4jProperties);
			
			// This Errorhandling should be uncommented for testing versions that should report errors:
			/*if (!LoggingSystem.alreadyInitiated()) {
				LoggingSystem.init(
						"/org/mutoss/gui/commons-logging.properties",
						//log4jProperties,
						true,
						false,
						new ApplicationLog());
				ErrorHandler.init("rohmeyer@small-projects.de", "http://www.algorithm-forge.com/report/bugreport.php", true, true, ErrorDialogChooseLevel.class);
			} */

			// Java 7 does not respect system property "sun.awt.exception.handler".
			// Eventually this fix should be included in afcommons.
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {		
					Thread.currentThread().setUncaughtExceptionHandler(new DefaultExceptionHandler());
				}
			});
		}
		
		if (Configuration.getInstance().getGeneralConfig().failSafeMode()) {
			JOptionPane.showMessageDialog(this, "The R enginge will now be initialized.", "Initializing R engine", JOptionPane.INFORMATION_MESSAGE);
		}
		RControl.getRControl(true);
		
		if (Configuration.getInstance().getGeneralConfig().failSafeMode()) {
			JOptionPane.showMessageDialog(this, "R version will be inspected", "R version", JOptionPane.INFORMATION_MESSAGE);
		}
		/* Get and save R and gMCP version numbers */
		try {		
			Configuration.getInstance().getGeneralConfig().setRVersionNumber(RControl.getR().eval("paste(R.version$major,R.version$minor,sep=\".\")").asRChar().getData()[0]);
			Configuration.getInstance().getGeneralConfig().setVersionNumber(RControl.getR().eval("Crossover:::CrossoverVersion()").asRChar().getData()[0]);
			this.setTitle("CROSS-OVER DESIGN SEARCH TOOL "+Configuration.getInstance().getGeneralConfig().getVersionNumber());
		} catch (Exception e) {
			// This is no vital information and will fail for e.g. R 2.8.0, so no error handling here...
		}
		
		if (!Configuration.getInstance().getGeneralConfig().failSafeMode() 
				|| JOptionPane.showConfirmDialog(this, "Set Icon and Look'n'Feel?", "Icon and Look and Feel", JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION) {
			setIconImage((new ImageIcon(getClass().getResource("/org/mutoss/gui/rjavaicon64.png"))).getImage());
			try {
				setLooknFeel();
			} catch (Exception e1) {
				JOptionPane.showMessageDialog(this, "Font size and Look'n'Feel could not be restored.", "Error restoring Look'n'Feel", JOptionPane.ERROR_MESSAGE);
			}
		}
		
		glassPane = new InfiniteProgressPanel(this, "Calculating");
	    setGlassPane(glassPane);
	    glassPane.addAbortListener(this);
		
		dac = new DesignArchiveControl();
	    
		makeContent();
		setJMenuBar(new MenuBarCrossover(this));
		pack();

		 // Place the frame in the middle of the screen with a border of inset = 50 pixel.
		int inset = 50;
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		setBounds(inset, inset,
				screenSize.width  - inset*2,
				screenSize.height - inset*2);

		/* 
		 * We want to check for unsaved changes and eventually quit the R console as well, 
		 * so we implement the WindowListener interface and let windowClosing() do the work.
		 */
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		addWindowListener(this);
		
		if (!Configuration.getInstance().getGeneralConfig().failSafeMode() 
				|| JOptionPane.showConfirmDialog(this, "Show GPL?", "GPL", JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION) {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {				
					if (!Configuration.getInstance().getGeneralConfig().haveReadGPL()) {
						new TellAboutOnlineUpate(null);
						Configuration.getInstance().getGeneralConfig().setTellAboutCheckOnline(true);
					}
					/*new Thread(new Runnable() {
					public void run() {
						Thread.currentThread().setUncaughtExceptionHandler(new DefaultExceptionHandler());
						VersionComparator.getOnlineVersion();
					}
				}).start();*/				
				}
			});
		}
		
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {				
				try {
				    Thread.sleep(1000);
				} catch(InterruptedException ex) {
				    Thread.currentThread().interrupt();
				}
				lock = false;
			}
		});
		
		setVisible(true);
		
		//JOptionPane.showMessageDialog(this, "THIS PACKAGE IS AN ALPHA VERSION!\nIt may give you good designs or it may make mistakes.\nPlease check for every design whether it fits your needs.", "ALPHA VERSION", JOptionPane.WARNING_MESSAGE);
	}

	/** Two methods to set and save the default values */
	private String getValue(String key, String def) { return conf.getProperty(key, def); }
	private void setValue(String key, String value) { conf.setProperty(key, value); }

	DesignSelectionPanel designPanel;
	DesignInputPanel designInputPanel;
	public AlgorithmPanel algorithmPanel;
	JTabbedPane tabbedPane;
	
	private void makeContent() {
		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.BOTH;
		c.gridx=0; c.gridy=0; c.gridwidth = 1; c.gridheight = 1; c.ipadx=5; c.ipady=5; c.weightx=1; c.weighty=0;
		
		getContentPane().setLayout(new GridBagLayout());		
		
		getContentPane().add(getInterface(), c);
		
		c.gridx++;
		
		getContentPane().add(getModelPanel(), c);
		
		tabbedPane = new JTabbedPane();
		
		designPanel = new DesignSelectionPanel(this);
		designInputPanel = new DesignInputPanel(this);
		designInputPanel.addActionListener(this);
		algorithmPanel = new AlgorithmPanel(this);
		tabbedPane.addTab("Catalogue", designPanel);
		tabbedPane.addTab("Algorithm Search", algorithmPanel);
		tabbedPane.addTab("Input own design", designInputPanel);
		tabbedPane.addChangeListener(this);
		
		c.gridy++; c.gridx=0; c.gridwidth=2;	 
		c.weighty=1;
		getContentPane().add(tabbedPane, c);	
		
		loadDefaults();
		newModelSelected();
		
		spinnerT.addChangeListener(this);
		spinnerP.addChangeListener(this);
		spinnerS1.addChangeListener(this);
		spinnerS2.addChangeListener(this);
		jCBmodel.addActionListener(this);
		
		stateChanged(null);
	}
	
	JSpinner spinnerT;
	JSpinner spinnerP;
	JSpinner spinnerS1;
	JSpinner spinnerS2;
	
	final int maxSp = 99;
	
	public JPanel getInterface() {
		JPanel panel = new JPanel();
		String cols = "5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu, fill:pref:grow";
        String rows = "5dlu, pref, 5dlu, pref, 5dlu, pref, 5dlu";
        
        panel.setLayout(new FormLayout(cols, rows));
        CellConstraints cc = new CellConstraints();
		
		int row = 2;
		
		spinnerT = new JSpinner(new SpinnerNumberModel(4, 2, maxSp, 1));    	
    	
    	
    	panel.add(new JLabel("Number of treatments:"), cc.xyw(2, row, 3));
        panel.add(spinnerT, cc.xy(6, row));
		
        row+=2;
		
    	spinnerP = new JSpinner(new SpinnerNumberModel(4, 2, maxSp, 1));    	
    	
    	
    	panel.add(new JLabel("Number of periods:"), cc.xyw(2, row, 3));
        panel.add(spinnerP, cc.xy(6, row));		
		
        row+=2;
		
    	spinnerS1 = new JSpinner(new SpinnerNumberModel(4, 2, maxSp, 1));    	
    	
		
    	spinnerS2 = new JSpinner(new SpinnerNumberModel(8, 2, maxSp, 1));    	
    	
    	
    	panel.add(new JLabel("Number of sequences:"), cc.xy(2, row));
    	panel.add(new JLabel("Min:"), cc.xy(4, row));
        panel.add(spinnerS1, cc.xy(6, row));
        panel.add(new JLabel("Max:"), cc.xy(8, row));
        panel.add(spinnerS2, cc.xy(10, row));
        
		return panel;
	}
	
	JComboBox jCBmodel;
	JLabel pLabel = new JLabel("Further model parameters:");
	JTextField jtfParam = new JTextField("1");
	
	public JPanel getModelPanel() {
		CellConstraints cc = new CellConstraints();
		JPanel modelPanel = new JPanel();
		String cols = "5dlu, pref, 5dlu, fill:min:grow, 5dlu";
		String rows = "5dlu, pref, 5dlu, pref, 5dlu";

		modelPanel.setLayout(new FormLayout(cols, rows));

		int row = 2;

		jCBmodel = new JComboBox(CrossoverGUI.models);
		

		modelPanel.add(new JLabel("Model"), cc.xy(2, row));
		modelPanel.add(jCBmodel, cc.xy(4, row));

		row+=2;		
        
		modelPanel.add(pLabel, cc.xy(2, row));
		modelPanel.add(jtfParam, cc.xy(4, row));
		pLabel.setEnabled(false);
		jtfParam.setEnabled(false);
		jtfParam.getDocument().addDocumentListener(this);
		return modelPanel;
	}
	
	ClassConfig ac = new ClassConfig(Configuration.getInstance(), CrossoverGUI.class);
	
	public void loadDefaults() {
		//JOptionPane.showMessageDialog(this, "Loading Settings");
		spinnerS1.getModel().setValue(ac.getIntProperty("s1", 4));
		spinnerS2.getModel().setValue(ac.getIntProperty("s2", 8));
		spinnerT.getModel().setValue(ac.getIntProperty("t", 4));
		spinnerP.getModel().setValue(ac.getIntProperty("p", 4));
		jCBmodel.setSelectedIndex(ac.getIntProperty("CVPattern", 0));
		jtfParam.setText(ac.getProperty("modelParam", "1"));
	}
    
	public void saveDefaults() {
		//JOptionPane.showMessageDialog(this, "Saving Settings");
		ac.setIntProperty("s1", Integer.parseInt(spinnerS1.getModel().getValue().toString()));
		ac.setIntProperty("s2", Integer.parseInt(spinnerS2.getModel().getValue().toString()));
		ac.setIntProperty("t", Integer.parseInt(spinnerT.getModel().getValue().toString()));
		ac.setIntProperty("p", Integer.parseInt(spinnerP.getModel().getValue().toString()));
		ac.setIntProperty("CVPattern", jCBmodel.getSelectedIndex());
		ac.setProperty("modelParam", jtfParam.getText());
	}
	
	/**
	 * This function handles events from buttons, check boxes and combo boxes. 
	 */
	public void actionPerformed(ActionEvent e) {
		saveDefaults(); //TODO: Should we call saveDefaults at other places (stateChanged - check, Model parameter Text field?)?
		if (e.getSource() == jCBmodel) {
			newModelSelected();
		}
	}
	
	public double getParam() {
		return Double.parseDouble(jtfParam.getText());
	}
	 
	private void newModelSelected() {
		algorithmPanel.createEffPanel();
		if (jCBmodel.getSelectedIndex()==PLACEBOMODEL) {
			jtfParam.setEnabled(true);
			pLabel.setEnabled(true);
			pLabel.setText("Number of placebo treatments:");
		} else if (jCBmodel.getSelectedIndex()==PROPORTIONALMODEL) {
			jtfParam.setEnabled(true);
			pLabel.setEnabled(true);
			pLabel.setText("Proportionality parameter:");
		} else {
			jtfParam.setEnabled(false);
			pLabel.setEnabled(false);
			pLabel.setText("Further model parameters:");
		}
		designPanel.valueChanged(null);	
		int i = designPanel.designTable.getSelectedRow();		
		designPanel.designTable.getModel().fireTableDataChanged();
		if (i != -1) {
			designPanel.designTable.setRowSelectionInterval(i, i);
		}
	}

	/**
	 * Shows a file from a package
	 * @param s String like 
	 */
	public void showFile(String s) {
		File f = new File(RControl.getR().eval("system.file(\""+s+"\", package=\"Crossover\")").asRChar().getData()[0]);
		if (!f.exists()) {
			throw new RuntimeException("This is strange. The file \""+s+"\" could not be found.");
		} else {
			try {	
				Method main = Class.forName("java.awt.Desktop").getDeclaredMethod("getDesktop");
				Object obj = main.invoke(new Object[0]);
				Method second = obj.getClass().getDeclaredMethod("open", new Class[] { File.class }); 
				second.invoke(obj, f);
			} catch (Exception exc) {			
				try {
					if (OSTools.isWindows()) {
						Process p;							
						p = Runtime.getRuntime().exec("rundll32 " +
								"url.dll,FileProtocolHandler " + f.getAbsolutePath());
						p.waitFor();
					} else {
						JOptionPane.showMessageDialog(this, "Please open and read the following file:\n"+f.getAbsolutePath(), "Could not find appropriate viewer", JOptionPane.WARNING_MESSAGE);
					}
				} catch (Exception e1) {					
					e1.printStackTrace();
					JOptionPane.showMessageDialog(this, "Please open and read the following file:\n"+f.getAbsolutePath(), "Could not find appropriate viewer", JOptionPane.WARNING_MESSAGE);
				}

			}
		}
	}

	/**
	 * WindowListener methods - the first one closes the R console if  
	 */
	public void windowClosing(WindowEvent e) {
		if (dac.designsChanged) {
			int answer = JOptionPane.showConfirmDialog(this, "There are unsaved designs. Should they be saved?", "Unsaved designs", JOptionPane.YES_NO_OPTION);
			if (answer==JOptionPane.YES_OPTION) {
				dac.save(this);
				// It is not save to continue if saving is still in progress.
				return;
			}
		}
		if (RControl.getR().eval("exists(\".isBundle\")").asRLogical().getData()[0]) {
			RControl.getR().eval("q(save=\"no\")");
		} else {
			if (System.getProperty("eclipse") != null) {
				System.exit(0);
			} else {
				dispose();
			}
		}
	}
	public void windowActivated(WindowEvent e) {}
	public void windowClosed(WindowEvent e) {}
	public void windowDeactivated(WindowEvent e) {}
	public void windowDeiconified(WindowEvent e) {}
	public void windowIconified(WindowEvent e) {}
	public void windowOpened(WindowEvent e) {}
	

	public void stateChanged(ChangeEvent e) {
		if (e!=null) saveDefaults();
		//if (e.getSource()==tabbedPane)		
		if (e!=null && (e.getSource()==spinnerS1 || e.getSource()==spinnerS2)) {
			int s1 = Integer.parseInt(spinnerS1.getModel().getValue().toString());
			int s2 = Integer.parseInt(spinnerS2.getModel().getValue().toString());
			int s = Integer.parseInt(algorithmPanel.spinnerS.getModel().getValue().toString());
			if (s1>s2) {
				if (e.getSource()==spinnerS1) {
					spinnerS2.getModel().setValue(s1);
				} else {
					spinnerS1.getModel().setValue(s2);
				}
			}
			if (s<s1) {
				algorithmPanel.spinnerS.getModel().setValue(s1);
			} else if (s2<s) {
				algorithmPanel.spinnerS.getModel().setValue(s2);
			}
		}
		if (tabbedPane.getSelectedIndex()!=0) return;
		glassPane.start();

		SwingWorker<Void, Void> worker = new FillTableWorker(this);
		worker.execute();		
	}

	/**
	 * Method to call from R
	 */
	public static void startGUI() {
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				new CrossoverGUI();
			}
		});		
	}	
	
	/**
	 * Just for me for the development
	 * @param args
	 */
	public static void main(String[] args) {
		RControl.getR().eval("Sys.setenv(\"JAVAGD_CLASS_NAME\"=\"org/mutoss/gui/JavaGD\")");
		new CrossoverGUI();
	}

	public void valueChanged(ListSelectionEvent e) {
		// TODO Auto-generated method stub
		
	}

	public void abort() {
		// TODO Auto-generated method stub
		
	}

	public String getParameters() {
		int t = Integer.parseInt(spinnerT.getModel().getValue().toString());
		int p = Integer.parseInt(spinnerP.getModel().getValue().toString());
		return "p="+p+", v="+t;
	}
	
	private void setLooknFeel() throws ClassNotFoundException, IllegalAccessException,
	InstantiationException, UnsupportedLookAndFeelException {
		UIManager.setLookAndFeel(Configuration.getInstance().getJavaConfig().getLooknFeel());
		WidgetFactory.setFontSizeGlobal(Configuration.getInstance().getGeneralConfig().getFontSize());
		SwingUtilities.updateComponentTreeUI(this);
	}

	private boolean lock = true;
	
	/**
	 * We must make sure that we return on the same thread that we entered from R.
	 * FillTableWorker must therefore wait until CrossoverGUI is finished.
	 * See for example: https://github.com/s-u/rJava/issues/18
	 * @return
	 */
	public synchronized boolean getLock() {
		if (!lock) {
			lock = true;
			return true;
		}
		return false;
	}
	
	public void resetLock() {
		lock = false;
	}

	public int getModel() {		
		return jCBmodel.getSelectedIndex()+1;
	}

	public void insertUpdate(DocumentEvent e) {
		try {
			newModelSelected();
		} catch(Exception ex) {
			
		}		
	}

	public void removeUpdate(DocumentEvent e) {
		try {
			newModelSelected();
		} catch(Exception ex) {
			
		}		
	}

	public void changedUpdate(DocumentEvent e) {
		try {
			newModelSelected();
		} catch(Exception ex) {
			
		}		
	}

}
